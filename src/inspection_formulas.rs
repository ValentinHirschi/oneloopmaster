//! Inspection-only, algebraically compact forms of the same sector formulas.
//!
//! The numerical maps retain their conditioned half-turn representation. For
//! inspection, carry the physical argument and its continued logarithm instead:
//! multiplication/division then add/subtract logarithms, without duplicating a
//! sign-normalized value and its parity throughout every downstream formula.
#![allow(clippy::duplicate_mod)] // Deliberately compile sectors against two sheet types.

use crate::{
    Atom, AtomCore, FunctionMap, LaurentSeries, MappedLaurentSeries, OneLoopExpressions, Symbol,
    any_nonzero, grouped, if_nonzero, if_nonzero_else, log_over_one_minus, native_function_options,
    off_real_axis, physical_log,
};

#[path = "inspection_sheet.rs"]
mod sheet_exact;
use sheet_exact::{SheetAtom, SheetAtom as ExactSheetAtom};
use sheet_exact::{
    divided_difference as sheet_dilog_divided_difference,
    divided_difference as exact_sheet_difference,
};

// Compile the sector formulas twice against the two sheet representations;
// there is one source of truth for kinematics, roots and sector dispatch.
#[path = "triangle.rs"]
#[allow(dead_code)] // The shared source also contains public numerical wrappers.
mod triangle;
use triangle::*;
#[path = "box_complex.rs"]
mod box_complex;
#[path = "box_integral.rs"]
#[allow(dead_code)] // Only sector construction is used in this inspection map.
mod box_integral;
#[path = "triangle_hv.rs"]
mod triangle_hv;

pub(crate) fn definitions() -> &'static FunctionMap {
    static MAP: std::sync::OnceLock<FunctionMap> = std::sync::OnceLock::new();
    MAP.get_or_init(|| {
        let mut map = FunctionMap::new();
        sheet_exact::register(&mut map);
        triangle_hv::register(&mut map);
        box_complex::register(&mut map);
        register_triangles(&mut map);
        box_integral::register_boxes(&mut map);
        crate::masters::register_masters(&mut map);
        map
    })
}

#[cfg(test)]
mod tests {
    pub(crate) use crate::tests::large_stack;
}
