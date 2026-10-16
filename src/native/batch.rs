//! Storage for emitted instruction-major kernels. No evaluator or IR is retained.

use super::primitives::{C, NativeFloat};
use symbolica::domains::float::FloatLike;

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
