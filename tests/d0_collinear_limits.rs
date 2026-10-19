//! D0 with at most two (adjacent) massive lines next to a collinear limit: an
//! invariant between two massless lines that is small but nonzero.
//!
//! The finite massless, one-mass and two-adjacent-mass formulas (OneLOop's
//! boxf0, boxf1 and boxf2) take dilogarithms of `X r34/r13` and `X r24/r12`
//! for the roots `X` of a quadratic, which tends to `(r34 X - r13)(r24 X -
//! r12)` as the invariant `r23` between the massless lines vanishes (or `r14`
//! without masses): one argument of each pair approaches one. On sheets with
//! windings the dilogarithm carries `i pi n log(1 - x)` there, the collinear
//! logarithm, and `1 - x` formed from a rounded argument lost `|1 - x|`
//! relative digits: the boxes below kept 1e-9 to 3e-7 at invariants of 1e-6
//! to 1e-10 and were O(1) wrong at 1e-14 (OneLOop's double precision loses the
//! same). Where the leading coefficient is the small one (`r34` or `r24`
//! small), the textbook small root lost them as well. The roots now follow
//! `solabc` (the larger directly, the other from Vieta's product), and `1 - x`
//! comes from the roots' offsets from `r13/r34` and `r12/r24`, where the
//! quadratic is the small invariant times a closed form.
//!
//! Next to one, sheets with two or more windings jump across the unit circle;
//! the side is decided from `1 - x` too. Next to `r13 r24 = r12 r34` as well,
//! the discriminant takes an exact form.
#![cfg(feature = "generated-evaluators")]
use oneloop::{
    EvaluationBackend, NativeEvaluator, PrecisionEvaluator, ScalarIntegral, evaluate_with_backend,
};
use symbolica::domains::float::{Complex, Float, RealLike};

type C = Complex<f64>;

/// Reference precision of the binary64 checks.
const BITS: u32 = 192;

fn on_stack(test: impl FnOnce() + Send + 'static) {
    std::thread::Builder::new()
        .stack_size(128 * 1024 * 1024)
        .spawn(test)
        .unwrap()
        .join()
        .unwrap();
}

/// `kind mu² p.. (re im)..` rows as evaluator inputs (mu² last).
fn input(row: &str) -> Vec<C> {
    let v: Vec<f64> = row.split_whitespace().map(|x| x.parse().unwrap()).collect();
    assert_eq!(v.len(), 16);
    let mut args: Vec<C> = v[2..8].iter().map(|&x| C::new(x, 0.0)).collect();
    args.extend(v[8..].chunks_exact(2).map(|m| C::new(m[0], m[1])));
    args.push(C::new(v[1], 0.0));
    args
}

fn high(args: &[C], bits: u32) -> Vec<Complex<Float>> {
    args.iter()
        .map(|z| Complex::new(Float::with_val(bits, z.re), Float::with_val(bits, z.im)))
        .collect()
}

fn relative(value: C, reference: C) -> f64 {
    let difference = value - reference;
    difference.re.hypot(difference.im) / reference.re.hypot(reference.im)
}

/// A one-mass box with a complex mass and p1² small, a one-mass box with a
/// real mass and p2² small, a two-adjacent-mass box with p4² small, and
/// massless boxes with p2² and p4² small, each at the reported invariant and
/// at ±10^-k for k = 2, 4, .., 14.
fn collinear_rows() -> Vec<String> {
    let boxes = [
        (
            "4 1.0 {} 16.245027303866983 -1.165189564810856 7.097305775236421 \
             -17.17229590533035 -16.317103653229925 0 0 0 0 0 0 1.0 -0.001",
            -6.657258425283804e-9,
        ),
        (
            "4 1.0 17.460951017893965 {} -18.651874759191873 4.55723164310599 \
             17.743781622935032 -17.07742653020065 0 0 0 0 0 0 1.0 0",
            -1.2687298815009212e-10,
        ),
        (
            "4 1.0 -2.754641313889903 11.792319055714035 6.738076671399869 {} \
             13.572798592328985 12.013797323318206 0 0 2.25 -0.00225 1.0 0 0 0",
            -4.460337301814886e-6,
        ),
        (
            "4 1.0 15.707984658659143 {} 4.9990570087914605 -5.6326442138090105 \
             -9.139526304659672 7.401553608499611 0 0 0 0 0 0 0 0",
            7.780146369002195e-14,
        ),
        (
            "4 1.0 2.7830216830853693 1.3231870761018882 -4.416237022072643 {} \
             5.690610472327094 8.050027823332083 0 0 0 0 0 0 0 0",
            -2.0329177941275214e-12,
        ),
    ];
    let mut rows = Vec::new();
    for (template, reported) in boxes {
        let row = |invariant: f64| template.replacen("{}", &format!("{invariant:e}"), 1);
        rows.push(row(reported));
        for k in (2..=14).step_by(2) {
            for sign in [-1.0, 1.0] {
                rows.push(row(sign * 10f64.powi(-k)));
            }
        }
    }
    rows
}

/// Rows whose binary64 routes are within their tolerance of native code at
/// 192 bits, or the failures.
fn binary64_failures(rows: impl IntoIterator<Item = (String, f64)>) -> Vec<String> {
    let mut precise =
        NativeEvaluator::<Float>::with_binary_precision(ScalarIntegral::D0, BITS).unwrap();
    let mut native = [true, false].map(|handwritten| {
        let mut evaluator = NativeEvaluator::<f64>::new(ScalarIntegral::D0).unwrap();
        evaluator.set_handwritten(handwritten);
        (handwritten, evaluator)
    });
    let mut failures = Vec::new();
    for (row, tolerance) in rows {
        let args = input(&row);
        let mut out =
            core::array::from_fn::<_, 3, _>(|_| Complex::new(Float::new(BITS), Float::new(BITS)));
        precise.evaluate(&high(&args, BITS), &mut out).unwrap();
        let reference = C::new(out[0].re.to_f64(), out[0].im.to_f64());
        let mut routes = Vec::new();
        for (handwritten, evaluator) in &mut native {
            let mut output = [C::new(f64::NAN, f64::NAN); 3];
            evaluator.evaluate(&args, &mut output).unwrap();
            routes.push((format!("native (hand-written: {handwritten})"), output[0]));
        }
        for backend in [EvaluationBackend::Expression, EvaluationBackend::SymJit] {
            let mut output = [C::new(f64::NAN, f64::NAN); 3];
            evaluate_with_backend(ScalarIntegral::D0, &args, &mut output, backend).unwrap();
            routes.push((backend.name().to_string(), output[0]));
        }
        for (route, value) in routes {
            let error = relative(value, reference);
            if !(error <= tolerance) {
                failures.push(format!("{error:.1e} {route} {row}: {value:?} vs {reference:?}"));
            }
        }
    }
    failures
}

/// Every binary64 route is within 1e-12 of native code at 192 bits (at most
/// 2e-13 here, SymJIT). Before, the reported boxes kept 2.9e-7, 4.2e-9 and
/// 1.1e-9, the scans up to 0.16 at 1e-14, and the massless boxes 4e-3 and
/// 1.4e-4 at their reported invariants.
#[test]
fn collinear_limits_in_binary64() {
    on_stack(|| {
        let failures = binary64_failures(collinear_rows().into_iter().map(|row| (row, 1e-12)));
        assert!(failures.is_empty(), "{failures:#?}");
    });
}

/// Collinear limits next to `r13 r24 = r12 r34`, where the quadratic's
/// discriminant is small and the dilogarithm terms grow as its inverse: a
/// one-mass box (p1² = 1.2e-7, kept 3.1e-5 before `b² - 4ac` took its exact
/// form, now within 1e-7), a two-adjacent-mass box with `r12 r34/(r13 r24) -
/// 1 = 2e-3` and p2² from 3.3 to 3.3e-14 (up to 1.8e-11 before, 9e-12 with
/// the exact discriminant alone, OneLOop's double precision 4e-13 at 3.3e-2
/// and 5e-3 at 3.3e-14; now within 1e-12, once both terms take one rounded
/// `r13 r24 - r12 r34`) and a massless box (p4²/mu² = -8.9e-12 with
/// `r12 r34/(r13 r24) - 1 = 2e-3`: 2.2e-2, now 4e-10 in native code and 5e-11
/// interpreted, the rounding of `r12 r34 - r13 r24` amplified; OneLOop's
/// double precision keeps 1.3e-2). With a massless box whose small p3² puts
/// the arguments far from one (the complements then carry no digits and must
/// not replace the values: 2.5e-12 when they did).
#[test]
fn next_to_a_vanishing_discriminant_in_binary64() {
    on_stack(|| {
        let rows = [
            (
                "4 1.0 1.2274203147474157e-07 40.19549357093282 -19.130300380643092 \
                 -1.733668902000268 13.95191697639357 -13.761223519470311 \
                 0 0 0 0 0 0 4.660339257746877 -0.004660339257746877",
                1e-7,
            ),
            (
                "4 0.10883312701972196 -0.18092389442831602 0.2304457119477068 \
                 2.1544187276504205 -9.654971344635188e-13 -1.497147830303522 \
                 0.2598589213555622 0 0 0 0 0 0 0 0",
                1e-9,
            ),
            (
                "4 19.959692369608153 -47.01132214535428 232.1725235092141 \
                 2.3532575220095417e-08 177.08313736745058 126.892513125016 \
                 320.89830093324014 0 0 0 0 0 0 0 0",
                1e-12,
            ),
        ];
        let two_masses = [0.0, 1.0, 2.0, 3.0, 3.5, 4.0, 5.0, 6.0, 8.0, 10.0, 12.0, 14.0].map(|k| {
            let p2 = 3.2740966864054275 * 10f64.powf(-k);
            let row = format!(
                "4 1.0 12.499100188334964 {p2:e} -17.603632557714292 2.901811859132522 \
                 15.904786899832231 -12.102786824580797 2.0499517118506425 0 0 0 0 0 \
                 4.8993952883173115 -0.2729525823783077"
            );
            (row, 1e-12)
        });
        let rows = rows.map(|(row, tolerance)| (row.to_string(), tolerance));
        let failures = binary64_failures(rows.into_iter().chain(two_masses));
        assert!(failures.is_empty(), "{failures:#?}");
    });
}

/// Random boxes (two adjacent masses, one real, with p2² = 7.6e-6 and
/// p4² = -9.8e-4; one complex mass with p1² = -9.3e-10) against OneLOop in
/// quadruple precision, in `Float` from 128 to 512 bits: native code and the
/// expression interpreter agree to 1e-15 at every precision. Native code was
/// O(1) wrong at 256 bits in all three and at 512 bits in two while the
/// unit-circle test signed a complex `1 - x conj(x)` (before 95f8107); the
/// complements' own test of the side must keep them right. The second box
/// (complex m2², real m3²) is one where OneLOop's two-mass box is wrong
/// (0.0338 + 0.1177i in quadruple precision, which these expressions
/// followed): its reference is Denner and Dittmaier's two-mass formula at 40
/// digits, which COLLIER's COLI and DD branches in quadruple precision
/// reproduce to 1e-19 (tests/d0_two_mass_width.rs).
#[test]
fn collinear_limits_at_every_precision() {
    on_stack(|| {
        let boxes = [
            (
                "4 1.0 13.339906506288962 7.645101522539689e-06 -18.067033948403353 \
                 1.295639379991652 15.748871904577292 -11.968624438049492 \
                 5.446210836482341 0 0 0 0 0 3.1100583835633797 -0.3560079210336333",
                C::new(-1.06460045004548420e-1, 7.00607484001557635e-3),
            ),
            (
                "4 1.0 -16.009634082627393 8.876629941904643 13.631357011185429 \
                 -0.0009822783830535126 3.455457549607029 -12.274987566569111 \
                 0 0 2.11121522382484 -0.0021112152238248402 2.7738362015154285 0 0 0",
                C::new(-1.77404138237007902e-1, 1.17659338364871247e-1),
            ),
            (
                "4 1.0 -9.305228475118043e-10 18.653938132067623 -16.47442840986097 \
                 14.781965947552756 -14.841260471245121 -19.288917101786765 \
                 0 0 0 0 0 0 0.3272646975008138 -0.0003272646975008138",
                C::new(8.91847915200563029e-2, 5.38406032947774360e0),
            ),
        ];
        let mut failures = Vec::new();
        for bits in [128, 160, 192, 224, 256, 384, 512] {
            let mut native =
                NativeEvaluator::<Float>::with_binary_precision(ScalarIntegral::D0, bits).unwrap();
            let mut expression = [128, 256].contains(&bits).then(|| {
                PrecisionEvaluator::with_binary_precision_and_backend(
                    ScalarIntegral::D0,
                    bits,
                    EvaluationBackend::Expression,
                )
                .unwrap()
            });
            for (row, expected) in boxes {
                let args = high(&input(row), bits);
                let mut routes = Vec::new();
                let fresh = || {
                    core::array::from_fn::<_, 3, _>(|_| {
                        Complex::new(Float::new(bits), Float::new(bits))
                    })
                };
                let mut output = fresh();
                native.evaluate(&args, &mut output).unwrap();
                routes.push(("native", output));
                if let Some(expression) = &mut expression {
                    let mut output = fresh();
                    expression.evaluate(&args, &mut output).unwrap();
                    routes.push(("expression", output));
                }
                for (route, output) in routes {
                    let value = C::new(output[0].re.to_f64(), output[0].im.to_f64());
                    let error = relative(value, expected);
                    if !(error <= 1e-15) {
                        failures.push(format!("{error:.1e} {route} {bits} bits {row}: {value:?}"));
                    }
                }
            }
        }
        assert!(failures.is_empty(), "{failures:#?}");
    });
}
