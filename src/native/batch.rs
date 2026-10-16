//! Storage for emitted instruction-major kernels. No evaluator or IR is retained.

use super::primitives::{C, NativeFloat};
use symbolica::domains::float::FloatLike;

/// Reusable slot-major storage for generated branch-free batch kernels.
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
}

impl<T: NativeFloat> BatchWorkspace<T> {
    pub fn new(prototype: T) -> Self {
        Self {
            prototype: C::new(prototype.zero(), prototype.zero()),
            values: Vec::new(),
            capacity: 0,
            slots: 0,
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

#[cfg(test)]
mod tests {
    use super::*;
    use symbolica::domains::float::Float;

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
