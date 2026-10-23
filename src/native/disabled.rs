//! Preserve the explicit native API so callers receive a useful configuration error.
use super::{C, NativeFloat};
use crate::ScalarIntegral;
use std::{convert::Infallible, marker::PhantomData};

/// The direct Rust evaluator requires the `generated-evaluators` Cargo feature.
/// Without it, constructors return an error; use the default evaluation API instead.
#[derive(Clone)]
pub struct NativeEvaluator<T: NativeFloat> {
    unavailable: Infallible,
    numeric_type: PhantomData<T>,
}

impl<T: NativeFloat> NativeEvaluator<T> {
    pub fn new(family: ScalarIntegral) -> Result<Self, String> {
        Self::with_binary_precision(family, T::DEFAULT_BITS)
    }

    pub fn with_binary_precision(_family: ScalarIntegral, _bits: u32) -> Result<Self, String> {
        Err("native backend requires the generated-evaluators Cargo feature; use the expression backend or rebuild with that feature".into())
    }

    pub fn family(&self) -> ScalarIntegral {
        match self.unavailable {}
    }

    pub fn binary_precision(&self) -> u32 {
        match self.unavailable {}
    }

    pub fn with_pole_selection(self, _selection: crate::PoleSelection) -> Self {
        match self.unavailable {}
    }

    pub fn set_pole_selection(&mut self, _selection: crate::PoleSelection) {
        match self.unavailable {}
    }

    pub fn pole_selection(&self) -> crate::PoleSelection {
        match self.unavailable {}
    }

    pub fn evaluate(&mut self, _input: &[C<T>], _output: &mut [C<T>]) -> Result<(), String> {
        match self.unavailable {}
    }

    pub fn evaluate_batch(
        &mut self,
        _input: &[C<T>],
        _output: &mut [C<T>],
        _rows: usize,
    ) -> Result<(), String> {
        match self.unavailable {}
    }
}
