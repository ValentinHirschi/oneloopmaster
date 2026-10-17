//! Storage for emitted instruction-major kernels. No evaluator or IR is retained.

use super::primitives::{C, NativeFloat};
use std::num::NonZeroUsize;
use symbolica::domains::float::FloatLike;

/// Explicit execution policy for an emitted row-major batch kernel.
///
/// The default retains one generated call for the complete batch. Fixed blocks
/// keep each row intact, including every flattened Dual component. This is an
/// execution/storage choice; it does not change the emitted arithmetic.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub enum BatchExecutionPolicy {
    #[default]
    Whole,
    FixedRows(NonZeroUsize),
}

impl BatchExecutionPolicy {
    /// Invoke an unchanged emitted batch entry point on complete row slices.
    /// `dimensions` contains the input and output component counts per row,
    /// including all derivative components when the graph is Dual-transformed.
    ///
    /// Both complete buffers and their products are validated before workspace
    /// mutation or any output write. Empty batches still make one empty call.
    /// Fixed blocks have natural tails, so a final single row uses the emitted
    /// scalar path. Each generated call initializes its own lazy branch masks.
    ///
    /// An oversized workspace is released when applying a fixed bound, including
    /// after prior whole-batch calls or a larger bound. Subsequent generated calls
    /// see at most the requested row count. The callback must use that row count
    /// with this workspace; directly calling `prepare` can bypass the policy.
    /// This bounds the retained row stride, not a kernel's number of slots.
    ///
    /// Only this policy's complete-buffer checks are transactional. A callback
    /// error is returned directly without rollback: the failing block can have
    /// written output, as can earlier blocks. Callers requiring whole-batch
    /// domain validation must perform it before this method.
    pub fn evaluate<T: NativeFloat>(
        self,
        input: &[C<T>],
        output: &mut [C<T>],
        rows: usize,
        dimensions: (usize, usize),
        workspace: &mut BatchWorkspace<T>,
        mut evaluate: impl FnMut(
            &[C<T>],
            &mut [C<T>],
            usize,
            &mut BatchWorkspace<T>,
        ) -> Result<(), String>,
    ) -> Result<(), String> {
        let (inputs, outputs) = dimensions;
        let input_len = rows
            .checked_mul(inputs)
            .ok_or("batch input size overflow")?;
        let output_len = rows
            .checked_mul(outputs)
            .ok_or("batch output size overflow")?;
        if input.len() != input_len || output.len() != output_len {
            return Err("batch buffer shape mismatch".into());
        }
        let Self::FixedRows(block) = self else {
            return evaluate(input, output, rows, workspace);
        };
        let block = block.get();
        workspace.release_oversized_rows(block);
        if rows == 0 {
            return evaluate(input, output, 0, workspace);
        }
        let mut first = 0;
        while first < rows {
            let count = (rows - first).min(block);
            let end = first + count;
            evaluate(
                &input[first * inputs..end * inputs],
                &mut output[first * outputs..end * outputs],
                count,
                workspace,
            )?;
            first = end;
        }
        Ok(())
    }
}

/// Borrowed numerical storage, parent stride, active lanes and control storage.
#[doc(hidden)]
pub type ControlledBatchBuffers<'a, T> = (&'a mut [C<T>], usize, &'a [usize], &'a mut [usize]);

/// Reusable slot-major storage for generated batch kernels.
///
/// Construct once with an explicit precision prototype, including for graphs
/// without inputs or constants. This prototype only initializes scratch cells;
/// input copies and computed values keep their own precision. Capacity is a
/// high-water mark and is retained for smaller batches and other kernels.
#[derive(Clone)]
pub struct BatchWorkspace<T: NativeFloat> {
    prototype: C<T>,
    values: Vec<C<T>>,
    capacity: usize,
    slots: usize,
    lanes: Vec<usize>,
    control: Vec<usize>,
    frames: usize,
}

impl<T: NativeFloat> BatchWorkspace<T> {
    pub fn new(prototype: T) -> Self {
        Self {
            prototype: C::new(prototype.zero(), prototype.zero()),
            values: Vec::new(),
            capacity: 0,
            slots: 0,
            lanes: Vec::new(),
            control: Vec::new(),
            frames: 0,
        }
    }

    /// Largest lane count reserved so far. Single-row kernels need no storage.
    pub fn capacity(&self) -> usize {
        self.capacity
    }

    /// Number of retained numerical cells (not the number of instructions).
    pub fn storage_len(&self) -> usize {
        self.values.len()
    }

    fn release_oversized_rows(&mut self, maximum: usize) {
        if self.capacity > maximum {
            self.values = Vec::new();
            self.lanes = Vec::new();
            self.control = Vec::new();
            self.capacity = 0;
            self.slots = 0;
            self.frames = 0;
        }
    }

    /// Reserve disjoint root/helper slots and borrow their storage and stride.
    /// Generated kernels overwrite each used cell before reading it. Changing
    /// stride therefore does not require repacking old numerical values.
    #[doc(hidden)]
    pub fn prepare(&mut self, rows: usize, slots: usize) -> Result<(&mut [C<T>], usize), String> {
        // Empty kernels/batches must not inflate the high-water stride: a
        // zero-slot call with usize::MAX rows has no storage requirements.
        if rows == 0 || slots == 0 {
            return Ok((&mut self.values, self.capacity));
        }
        let capacity = self.capacity.max(rows);
        let slots = self.slots.max(slots);
        let length = capacity
            .checked_mul(slots)
            .ok_or("batch workspace size overflow")?;
        // Vec's byte-size limit can be exceeded even without element overflow.
        if length > isize::MAX as usize / std::mem::size_of::<C<T>>().max(1) {
            return Err("batch workspace byte size overflow".into());
        }
        if length > self.values.len() {
            self.values
                .try_reserve(length - self.values.len())
                .map_err(|error| error.to_string())?;
            self.values.resize(length, self.prototype.clone());
        }
        self.capacity = capacity;
        self.slots = slots;
        Ok((&mut self.values, capacity))
    }

    /// Reserve numerical cells and disjoint active/resume columns for lazy
    /// generated functions. Frame indices follow the exported helper DAG.
    /// All dimensions are checked before mutating storage or its high-water marks.
    #[doc(hidden)]
    pub fn prepare_controlled(
        &mut self,
        rows: usize,
        slots: usize,
        frames: usize,
    ) -> Result<ControlledBatchBuffers<'_, T>, String> {
        if rows == 0 || slots == 0 {
            return Ok((&mut self.values, self.capacity, &[], &mut self.control));
        }
        let capacity = self.capacity.max(rows);
        let slots = self.slots.max(slots);
        let frames = self.frames.max(frames);
        let length = capacity
            .checked_mul(slots)
            .ok_or("batch workspace size overflow")?;
        let controls = frames
            .checked_mul(2)
            .and_then(|n| n.checked_mul(capacity))
            .ok_or("batch control workspace size overflow")?;
        if length > isize::MAX as usize / std::mem::size_of::<C<T>>().max(1)
            || capacity > isize::MAX as usize / std::mem::size_of::<usize>()
            || controls > isize::MAX as usize / std::mem::size_of::<usize>()
        {
            return Err("batch workspace byte size overflow".into());
        }
        // Reserve all arenas before changing their contents or metadata.
        self.values
            .try_reserve(length.saturating_sub(self.values.len()))
            .map_err(|error| error.to_string())?;
        self.control
            .try_reserve(controls.saturating_sub(self.control.len()))
            .map_err(|error| error.to_string())?;
        self.lanes
            .try_reserve(capacity.saturating_sub(self.lanes.len()))
            .map_err(|error| error.to_string())?;
        if length > self.values.len() {
            self.values.resize(length, self.prototype.clone());
        }
        self.control.resize(controls, 0);
        self.lanes.extend(self.lanes.len()..capacity);
        self.capacity = capacity;
        self.slots = slots;
        self.frames = frames;
        Ok((
            &mut self.values,
            capacity,
            &self.lanes[..rows],
            &mut self.control,
        ))
    }
}

/// Copy a value without fixing or increasing its tracked precision.
#[doc(hidden)]
#[inline(always)]
pub fn copy_value<T: NativeFloat>(target: &mut C<T>, source: &C<T>) {
    target.set_from(source);
}

/// Copy disjoint slot columns without cloning numerical payloads into a
/// temporary. Identical columns are a no-op; the checked layout is emitted.
#[doc(hidden)]
pub fn copy_column<T: NativeFloat>(
    values: &mut [C<T>],
    stride: usize,
    rows: usize,
    target: usize,
    source: usize,
) {
    if target == source || rows == 0 {
        return;
    }
    let target = target * stride;
    let source = source * stride;
    let (target, source) = if target > source {
        let (left, right) = values.split_at_mut(target);
        (&mut right[..rows], &left[source..source + rows])
    } else {
        let (left, right) = values.split_at_mut(source);
        (&mut left[target..target + rows], &right[..rows])
    };
    for (target, source) in target.iter_mut().zip(source) {
        copy_value(target, source);
    }
}

/// Copy only a caller's active lanes, retaining the parent numerical stride.
#[doc(hidden)]
#[inline]
pub fn copy_active_column<T: NativeFloat>(
    values: &mut [C<T>],
    stride: usize,
    lanes: &[usize],
    target: usize,
    source: usize,
) {
    if target == source || lanes.is_empty() {
        return;
    }
    let target = target * stride;
    let source = source * stride;
    let (target, source) = if target > source {
        let (left, right) = values.split_at_mut(target);
        (&mut right[..stride], &left[source..source + stride])
    } else {
        let (left, right) = values.split_at_mut(source);
        (&mut left[target..target + stride], &right[..stride])
    };
    for &lane in lanes {
        copy_value(&mut target[lane], &source[lane]);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use symbolica::domains::float::Float;

    #[test]
    fn fixed_policy_checks_complete_shapes_before_reset_or_writes() {
        let mut workspace = BatchWorkspace::new(0.0);
        workspace.prepare_controlled(128, 4, 3).unwrap();
        let before = (
            workspace.capacity(),
            workspace.values.clone(),
            workspace.control.clone(),
            workspace.lanes.clone(),
        );
        let policy = BatchExecutionPolicy::FixedRows(NonZeroUsize::new(32).unwrap());
        let input = vec![C::new(1.0, -0.0); 130];
        let mut output = vec![C::new(99.0, 7.0); 65];
        for (rows, dimensions, input_length) in [
            (65, (2, 1), 129),
            (64, (2, 1), 128),
            (usize::MAX, (2, 1), 0),
            (usize::MAX, (0, 2), 0),
        ] {
            assert!(
                policy
                    .evaluate(
                        &input[..input_length],
                        &mut output,
                        rows,
                        dimensions,
                        &mut workspace,
                        |_, _, _, _| panic!("invalid complete buffers must not invoke a block")
                    )
                    .is_err()
            );
            assert_eq!(output, vec![C::new(99.0, 7.0); 65]);
            assert_eq!(
                (
                    workspace.capacity(),
                    workspace.values.clone(),
                    workspace.control.clone(),
                    workspace.lanes.clone()
                ),
                before
            );
        }
    }

    #[test]
    fn policy_calls_whole_once_and_returns_block_errors_in_row_order() {
        let input: Vec<_> = (0..17).map(|i| C::new(i as f64, -0.0)).collect();
        for policy in [
            BatchExecutionPolicy::Whole,
            BatchExecutionPolicy::FixedRows(NonZeroUsize::new(8).unwrap()),
        ] {
            let mut workspace = BatchWorkspace::new(0.0);
            let mut output = vec![C::new(-99.0, 1.0); input.len()];
            let mut calls = Vec::new();
            let failing_call = if policy == BatchExecutionPolicy::Whole {
                1
            } else {
                2
            };
            let error = policy
                .evaluate(
                    &input,
                    &mut output,
                    17,
                    (1, 1),
                    &mut workspace,
                    |input, output, rows, _| {
                        calls.push((input[0].re as usize, rows));
                        if calls.len() == failing_call {
                            output[0] = input[0];
                            return Err("callback failure is preserved".into());
                        }
                        output.clone_from_slice(input);
                        Ok(())
                    },
                )
                .unwrap_err();
            assert_eq!(error, "callback failure is preserved");
            if policy == BatchExecutionPolicy::Whole {
                assert_eq!(calls, [(0, 17)]);
                assert_eq!(output[0], input[0]);
                assert!(output[1..].iter().all(|v| *v == C::new(-99.0, 1.0)));
            } else {
                assert_eq!(calls, [(0, 8), (8, 8)]);
                assert_eq!(&output[..9], &input[..9]);
                assert!(output[9..].iter().all(|v| *v == C::new(-99.0, 1.0)));
            }
        }
    }

    #[test]
    fn fixed_policy_preserves_rows_order_tails_and_releases_all_oversized_arenas() {
        let mut workspace = BatchWorkspace::new(0.0);
        for (policy, rows, expected_calls) in [
            (BatchExecutionPolicy::Whole, 129, vec![129]),
            (
                BatchExecutionPolicy::FixedRows(NonZeroUsize::new(32).unwrap()),
                65,
                vec![32, 32, 1],
            ),
            (
                BatchExecutionPolicy::FixedRows(NonZeroUsize::new(8).unwrap()),
                17,
                vec![8, 8, 1],
            ),
            (
                BatchExecutionPolicy::FixedRows(NonZeroUsize::new(8).unwrap()),
                3,
                vec![3],
            ),
            (
                BatchExecutionPolicy::FixedRows(NonZeroUsize::new(8).unwrap()),
                0,
                vec![0],
            ),
            (
                BatchExecutionPolicy::FixedRows(NonZeroUsize::new(16).unwrap()),
                33,
                vec![16, 16, 1],
            ),
            (BatchExecutionPolicy::Whole, 129, vec![129]),
            (
                BatchExecutionPolicy::FixedRows(NonZeroUsize::new(4).unwrap()),
                0,
                vec![0],
            ),
            (
                BatchExecutionPolicy::FixedRows(NonZeroUsize::new(4).unwrap()),
                1,
                vec![1],
            ),
        ] {
            // Five components per row model a flattened higher/mixed Dual graph.
            let input: Vec<_> = (0..rows * 5).map(|i| C::new(i as f64, -0.0)).collect();
            let mut output = vec![C::new(-99.0, 99.0); input.len()];
            let mut calls = Vec::new();
            policy
                .evaluate(
                    &input,
                    &mut output,
                    rows,
                    (5, 5),
                    &mut workspace,
                    |input, output, rows, workspace| {
                        calls.push(rows);
                        if rows > 1 {
                            let (_, _, lanes, _) = workspace.prepare_controlled(rows, 7, 2)?;
                            assert_eq!(lanes, (0..rows).collect::<Vec<_>>());
                        }
                        output.clone_from_slice(input);
                        Ok(())
                    },
                )
                .unwrap();
            assert_eq!(calls, expected_calls);
            for (actual, expected) in output.iter().zip(&input) {
                assert_eq!(actual.re.to_bits(), expected.re.to_bits());
                assert_eq!(actual.im.to_bits(), expected.im.to_bits());
            }
            if let BatchExecutionPolicy::FixedRows(bound) = policy {
                assert!(workspace.capacity <= bound.get());
                assert!(workspace.values.len() <= bound.get() * 7);
                assert!(workspace.lanes.len() <= bound.get());
                assert!(workspace.control.len() <= bound.get() * 4);
                if rows == 0 && bound.get() == 4 {
                    assert_eq!(workspace.values.capacity(), 0);
                    assert_eq!(workspace.lanes.capacity(), 0);
                    assert_eq!(workspace.control.capacity(), 0);
                }
            }
        }
    }

    #[test]
    fn fixed_policy_retains_prototype_and_input_precision() {
        let mut workspace = BatchWorkspace::new(Float::new(177));
        workspace.prepare_controlled(128, 9, 7).unwrap();
        let input: Vec<_> = (0..33)
            .map(|i| C::new(Float::with_val(53 + i, i), Float::with_val(256 - i, -1)))
            .collect();
        let mut output = vec![C::new(Float::new(96), Float::new(96)); input.len()];
        BatchExecutionPolicy::FixedRows(NonZeroUsize::new(8).unwrap())
            .evaluate(&input, &mut output, 33, (1, 1), &mut workspace, |input, output, rows, workspace| {
                if rows > 1 {
                    let (values, _, _, _) = workspace.prepare_controlled(rows, 2, 1)?;
                    assert!(values.iter().all(|v| v.re.get_precision() == 177 && v.im.get_precision() == 177));
                }
                for (target, source) in output.iter_mut().zip(input) {
                    copy_value(target, source);
                }
                Ok(())
            }).unwrap();
        for (actual, expected) in output.iter().zip(&input) {
            assert_eq!(actual, expected);
            assert_eq!(actual.re.get_precision(), expected.re.get_precision());
            assert_eq!(actual.im.get_precision(), expected.im.get_precision());
        }
        assert_eq!(workspace.capacity(), 8);
        assert_eq!(workspace.values.len(), 16);
        assert_eq!(workspace.lanes.len(), 8);
        assert_eq!(workspace.control.len(), 16);
        assert_eq!(workspace.frames, 1);
    }

    #[test]
    fn controlled_workspace_reuses_parent_stride_and_keeps_inactive_values() {
        let mut workspace = BatchWorkspace::new(Float::new(192));
        {
            let (values, stride, lanes, control) = workspace.prepare_controlled(8, 3, 2).unwrap();
            assert_eq!(lanes, &[0, 1, 2, 3, 4, 5, 6, 7]);
            assert_eq!(control.len(), 32);
            for lane in 0..stride {
                values[lane] = C::new(
                    Float::with_val(53 + lane as u32, lane),
                    Float::with_val(256, 1),
                );
                values[stride + lane] = C::new(Float::with_val(192, -99), Float::with_val(192, 0));
            }
            copy_active_column(values, stride, &[1, 5, 7], 1, 0);
            for lane in 0..stride {
                if [1, 5, 7].contains(&lane) {
                    assert_eq!(values[stride + lane], values[lane]);
                    assert_eq!(values[stride + lane].re.get_precision(), 53 + lane as u32);
                    assert_eq!(values[stride + lane].im.get_precision(), 256);
                } else {
                    assert_eq!(values[stride + lane].re, Float::with_val(192, -99));
                }
            }
        }
        let before = workspace.values.clone();
        let control = workspace.control.clone();
        for (rows, slots, frames) in [(usize::MAX, 3, 2), (8, 3, usize::MAX), (8, usize::MAX, 2)] {
            assert!(workspace.prepare_controlled(rows, slots, frames).is_err());
            assert_eq!(workspace.capacity(), 8);
            assert_eq!(workspace.values, before);
            assert_eq!(workspace.control, control);
        }
        workspace
            .prepare_controlled(usize::MAX, 0, usize::MAX)
            .unwrap();
        workspace
            .prepare_controlled(0, usize::MAX, usize::MAX)
            .unwrap();
        assert_eq!(workspace.capacity(), 8);
        workspace.prepare(16, 2).unwrap();
        let (_, stride, lanes, control) = workspace.prepare_controlled(3, 2, 1).unwrap();
        assert_eq!(stride, 16);
        assert_eq!(lanes, &[0, 1, 2]);
        assert_eq!(control.len(), 64);
        let mut fresh = BatchWorkspace::new(0.0);
        fresh.prepare_controlled(usize::MAX, 0, usize::MAX).unwrap();
        fresh.prepare_controlled(2, 1, 1).unwrap();
        assert_eq!(fresh.capacity(), 2);
    }

    #[test]
    fn workspace_dimensions_overflow_without_mutation_and_copies_preserve_precision() {
        let mut workspace = BatchWorkspace::new(Float::new(192));
        let too_many = isize::MAX as usize / std::mem::size_of::<C<Float>>() + 1;
        assert!(workspace.prepare(too_many, 1).is_err());
        assert_eq!(workspace.capacity(), 0);
        assert_eq!(workspace.storage_len(), 0);

        let (values, stride) = workspace.prepare(3, 4).unwrap();
        assert_eq!(stride, 3);
        values[0] = C::new(Float::with_val(53, 2), Float::with_val(320, 3));
        copy_column(values, stride, 1, 3, 0);
        assert_eq!(values[9], values[0]);
        assert_eq!(values[9].re.get_precision(), 53);
        assert_eq!(values[9].im.get_precision(), 320);
        copy_column(values, stride, 1, 0, 3);
        copy_column(values, stride, 1, 0, 0);
        workspace.prepare(usize::MAX, 0).unwrap();
        workspace.prepare(0, usize::MAX).unwrap();
        assert_eq!(workspace.capacity(), 3);
        assert_eq!(workspace.storage_len(), 12);
        let before = workspace.values.clone();
        for (rows, slots) in [(usize::MAX, 4), (3, usize::MAX), (isize::MAX as usize, 1)] {
            assert!(workspace.prepare(rows, slots).is_err());
            assert_eq!(workspace.values, before);
            assert_eq!(workspace.capacity(), 3);
            assert_eq!(workspace.storage_len(), 12);
        }
        workspace.prepare(1, 2).unwrap();
        assert_eq!(workspace.capacity(), 3);
        assert_eq!(workspace.storage_len(), 12);
        workspace.prepare(7, 2).unwrap();
        assert_eq!(workspace.capacity(), 7);
        assert_eq!(workspace.storage_len(), 28);
    }
}
