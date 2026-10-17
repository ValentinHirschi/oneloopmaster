//! Reusable definitions and evaluator construction for mixed integral expressions.
use super::*;
use std::sync::OnceLock;
use symbolica::{
    domains::{dual::DualNumberStructure, float::Complex, rational::Rational},
    evaluate::{Dualizer, EvaluationError, ExpressionEvaluator},
};

/// Native definitions shared by any number of triangle and box expressions.
/// Construct this once when building an amplitude, then pass its function map
/// to Symbolica or use [Self::evaluator] to compile a combination of integrals.
#[derive(Clone, Debug)]
pub struct OneLoopExpressions {
    map: &'static FunctionMap,
}

impl Default for OneLoopExpressions {
    fn default() -> Self {
        Self::new()
    }
}

impl OneLoopExpressions {
    /// Builds the native definitions once per process and reuses them thereafter.
    pub fn new() -> Self {
        Self {
            map: shared_definitions(),
        }
    }

    /// Byte size of each registered definition body, largest first. Building
    /// evaluators traverses these bodies, so their size bounds the time and
    /// memory of every evaluator construction.
    #[doc(hidden)]
    pub fn definition_sizes(&self) -> Vec<(String, usize)> {
        self.map.body_sizes()
    }

    /// Returns the shared transparent Symbolica definitions.
    pub fn function_map(&self) -> &symbolica::evaluate::FunctionMap {
        self.map.as_symbolica()
    }
    /// Transfers the definitions to a Symbolica evaluator builder.
    pub fn into_function_map(self) -> symbolica::evaluate::FunctionMap {
        self.map.as_symbolica().clone()
    }

    pub(crate) fn into_definitions(self) -> &'static FunctionMap {
        self.map
    }

    /// Constructs a triangle using this context's definitions.
    /// The same unsupported degenerate limits as [crate::c0] apply.
    pub fn c0(&self, p: [&Atom; 3], m: [&Atom; 3], mu: &Atom) -> LaurentSeries {
        triangle_series(p, m, mu)
    }

    /// Constructs a box using this context's definitions.
    /// The same coverage limitations as [crate::d0] apply.
    pub fn d0(&self, p: [&Atom; 6], m: [&Atom; 4], mu: &Atom) -> LaurentSeries {
        box_series(p, m, mu)
    }

    /// Compiles combinations of integral coefficients. The returned evaluator
    /// retains exact rational constants. See the README for stack requirements.
    pub fn evaluator(
        &self,
        expressions: &[Atom],
        parameters: &[Atom],
    ) -> Result<ExpressionEvaluator<Complex<Rational>>, EvaluationError> {
        compile(expressions, parameters, self.map, false)
    }

    /// Compiles flattened Dual components after preparing removable helper neighborhoods.
    ///
    /// The maximum total derivative order comes from `dual`, so the generated
    /// local polynomials cover every component of its shape. Inputs and outputs
    /// use Symbolica's ordinary argument-major/component-minor Dual layout.
    /// This currently preserves equal nonzero sheet quotient jets and the sheet
    /// `log(x)/(1-x)` limit at one. An unresolved `0/0` sheet quotient stays
    /// nonfinite instead of inheriting the historical scalar convention. Coincident
    /// mapped sheet dilogarithm differences retain independent argument jets,
    /// including the joint unit limit on the unwound sheet. At coincidence the
    /// sheet tags must match; zero arguments and wound unit points remain
    /// nonanalytic. Other triangle/box limits still need their own neighborhoods;
    /// this constructor does not establish a C0/D0 derivative domain. Callers
    /// must validate their domain independently.
    /// Shared-body Dual evaluation requires Symbolica's recursive vectorizer
    /// and control-flow scope fix for branches with multiple result components.
    pub fn dual_evaluator<D: DualNumberStructure>(
        &self,
        expressions: &[Atom],
        parameters: &[Atom],
        dual: D,
    ) -> Result<ExpressionEvaluator<Complex<Rational>>, String> {
        compile_dual(expressions, parameters, self.map, dual)
    }

    /// Compile a manually assembled amplitude with SymJIT O2. The returned
    /// evaluator supports scalar and row-major batched complex-f64 inputs.
    /// The exact/arbitrary-precision [Self::evaluator] API remains unchanged.
    pub fn jit_evaluator(
        &self,
        expressions: &[Atom],
        parameters: &[Atom],
    ) -> Result<JitEvaluator, String> {
        let exact = compile(expressions, parameters, self.map, true).map_err(|e| e.to_string())?;
        JitEvaluator::compile(&exact, parameters.len(), expressions.len())
    }
}

pub(crate) fn shared_definitions() -> &'static FunctionMap {
    initialization::ensure_symbolica_state();
    static DEFINITIONS: OnceLock<FunctionMap> = OnceLock::new();
    DEFINITIONS.get_or_init(|| {
        let mut map = FunctionMap::new();
        sheet_exact::register(&mut map);
        triangle_hv::register(&mut map);
        box_complex::register(&mut map);
        register_triangles(&mut map);
        register_boxes(&mut map);
        register_masters(&mut map);
        map
    })
}

// Pure expression construction: safe to call while the shared definition map
// is being initialized. Do not invoke OneLoopExpressions::new from here.
pub(crate) fn triangle_series(p: [&Atom; 3], m: [&Atom; 3], mu: &Atom) -> LaurentSeries {
    let args = p
        .into_iter()
        .chain(m)
        .chain([mu])
        .cloned()
        .collect::<Vec<_>>();
    let branches = core::array::from_fn(|sector| series_call("c0", sector, &args));
    with_triangle_normalization(select_three_series(m, branches))
}

pub(crate) fn box_series(p: [&Atom; 6], m: [&Atom; 4], mu: &Atom) -> LaurentSeries {
    // Original d0 changes the channel pair before mass-sector dispatch when
    // s12 or s23 vanishes. The finite sector formulas can divide by those
    // channels, so skipping this permutation can silently return zero.
    // Only the exact-zero case is selected here, not the reference's optional
    // tolerance-based "hard kinematics" warning/coercion.
    let abs = p.map(|v| Symbol::ABS.call((v,)));
    let min = |a: &Atom, b: &Atom| {
        if_nonzero_else(&sheet_exact::negative(&(a - b)), a.clone(), b.clone())
    };
    let min13 = min(&abs[0], &abs[2]);
    let min24 = min(&abs[1], &abs[3]);
    let min56 = min(&abs[4], &abs[5]);
    let choose13 =
        sheet_exact::negative(&(&min24 - &min13)) * sheet_exact::negative(&(&min56 - &min13));
    let choose24 =
        sheet_exact::negative(&(&min13 - &min24)) * sheet_exact::negative(&(&min56 - &min24));
    let hard = if_nonzero(p[4], if_nonzero(p[5], Atom::num(1)));
    let select = |ordinary: &Atom, paired13: &Atom, paired24: &Atom| {
        if_nonzero_else(
            &hard,
            ordinary.clone(),
            if_nonzero_else(
                &choose13,
                paired13.clone(),
                if_nonzero_else(&choose24, paired24.clone(), ordinary.clone()),
            ),
        )
    };
    let p13 = [4, 1, 5, 3, 0, 2];
    let p24 = [0, 5, 2, 4, 3, 1];
    let m13 = [0, 2, 1, 3];
    let m24 = [0, 1, 3, 2];
    let pp = core::array::from_fn::<_, 6, _>(|i| select(p[i], p[p13[i]], p[p24[i]]));
    let mm = core::array::from_fn::<_, 4, _>(|i| select(m[i], m[m13[i]], m[m24[i]]));
    let args = pp
        .iter()
        .chain(&mm)
        .chain([mu])
        .cloned()
        .collect::<Vec<_>>();
    select_four_series(
        mm.each_ref(),
        core::array::from_fn(|sector| series_call("d0", sector, &args)),
    )
}

fn series_call(family: &str, sector: usize, args: &[Atom]) -> LaurentSeries {
    let call = |coefficient| {
        symbolica::symbol!(format!("__olo_{family}_sector_{sector}_{coefficient}")).call(args)
    };
    LaurentSeries::new(call(0), call(1), call(2))
}

fn compile(
    expressions: &[Atom],
    parameters: &[Atom],
    map: &FunctionMap,
    jit: bool,
) -> Result<ExpressionEvaluator<Complex<Rational>>, EvaluationError> {
    let (expressions, parameters) = map.normalize_inputs(expressions, parameters, jit);
    Atom::evaluator_multiple(&expressions, &parameters)
        .function_map(
            if jit {
                map.as_jit_symbolica()
            } else {
                map.as_symbolica()
            }
            .clone(),
        )
        .direct_translation(true)
        // Keep the formulas' arithmetic order and skip whole-map Horner expansion.
        .horner_iterations(0)
        .build()
}

fn compile_dual<D: DualNumberStructure>(
    expressions: &[Atom],
    parameters: &[Atom],
    map: &FunctionMap,
    dual: D,
) -> Result<ExpressionEvaluator<Complex<Rational>>, String> {
    let order = dual
        .get_shape()
        .iter()
        .map(|component| {
            component.iter().try_fold(0usize, |sum, &degree| {
                sum.checked_add(degree).ok_or("Dual total order overflow")
            })
        })
        .collect::<Result<Vec<_>, _>>()?
        .into_iter()
        .max()
        .ok_or("Dual shape is empty")?;
    let map = sheet_exact::with_neighborhoods(map, order)?;
    compile(expressions, parameters, &map, false)
        .map_err(|error| error.to_string())?
        .vectorize(&Dualizer::new(dual, Vec::new()))
}

impl MappedLaurentSeries {
    /// Compiles the three coefficients using the included definitions.
    pub fn evaluator(
        &self,
        parameters: &[Atom],
    ) -> Result<ExpressionEvaluator<Complex<Rational>>, EvaluationError> {
        compile(self.coefficients(), parameters, self.function_map, false)
    }

    /// Compiles shape-derived helper neighborhoods and flattened Dual components.
    /// See [OneLoopExpressions::dual_evaluator] for the currently covered helper
    /// and the remaining triangle/box derivative-domain limitations.
    pub fn dual_evaluator<D: DualNumberStructure>(
        &self,
        parameters: &[Atom],
        dual: D,
    ) -> Result<ExpressionEvaluator<Complex<Rational>>, String> {
        compile_dual(self.coefficients(), parameters, self.function_map, dual)
    }

    /// Compile all three coefficients for scalar or batched SymJIT evaluation.
    pub fn jit_evaluator(&self, parameters: &[Atom]) -> Result<JitEvaluator, String> {
        let exact = compile(self.coefficients(), parameters, self.function_map, true)
            .map_err(|e| e.to_string())?;
        JitEvaluator::compile(&exact, parameters.len(), 3)
    }
}
