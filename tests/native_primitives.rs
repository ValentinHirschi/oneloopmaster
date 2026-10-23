//! Numerical primitives, including Symbolica's registered polylog callback.

#[path = "../src/native/primitives.rs"]
mod primitives;

use primitives::{C, NativeFloat};
use symbolica::domains::float::{DoubleFloat, Float, FloatLike, Real, RealLike, SingleFloat};

fn dd_float(value: DoubleFloat, bits: u32) -> Float {
    let pair = value.into_inner();
    Float::with_val(bits, pair.hi()) + Float::with_val(bits, pair.lo())
}

fn relative(actual: &Float, expected: &Float, bits: u32) {
    let error = (actual.clone() - expected).norm();
    let scale = if expected.is_fully_zero() {
        expected.one()
    } else {
        expected.norm()
    };
    let bound = scale * Float::with_val(bits, 0.5).pow(u64::from(bits - 16));
    assert!(
        error <= bound,
        "actual={actual:e}, expected={expected:e}, error={error:e}"
    );
}

fn axes<T: NativeFloat>() {
    let zero = T::zero_at(128);
    assert_eq!(zero.at_precision(128), zero);
    for sign in [-1, 1] {
        let signed_zero = if sign < 0 {
            -zero.clone()
        } else {
            zero.clone()
        };
        let negative = C::new(zero.from_i64(-4), signed_zero.clone());
        let root = primitives::sqrt(&negative);
        assert!(root.re.is_fully_zero());
        assert_eq!(root.im, zero.from_i64(2 * sign));
        let logarithm = primitives::log(&negative);
        assert_eq!(logarithm.re, zero.from_i64(4).log());
        assert_eq!(logarithm.im, zero.pi() * zero.from_i64(sign));
        let inverse = primitives::powi(&negative, -1, true);
        assert_eq!(inverse.re, -zero.one() / zero.from_i64(4));
        assert!(inverse.im.is_fully_zero());
        let half = C::new(zero.one() / zero.from_i64(2), zero.clone());
        assert_eq!(primitives::powf(&negative, &half, true), root);
        let minus_half = C::new(-half.re, zero.clone());
        let inverse_root = primitives::powf(&negative, &minus_half, true);
        let unit = primitives::mul(&root, &inverse_root);
        assert_eq!(unit.re, zero.one());
        assert!(unit.im.is_fully_zero());
        let positive = C::new(zero.from_i64(4), signed_zero.clone());
        let root = primitives::sqrt(&positive);
        assert_eq!(root.re, zero.from_i64(2));
        assert_eq!(root.im.negative_sign(), sign < 0);
        assert_eq!(primitives::log(&positive).im.negative_sign(), sign < 0);
        assert!(!primitives::truth(&C::new(
            signed_zero.clone(),
            signed_zero
        )));
    }
    let value = C::new(zero.from_i64(3), zero.from_i64(-4));
    let norm_squared = primitives::mul(&value, &primitives::conj(&value));
    assert_eq!(norm_squared.re, zero.from_i64(25));
    assert!(norm_squared.im.is_fully_zero());
    assert_eq!(
        primitives::abs(&value),
        C::new(zero.from_i64(5), zero.clone())
    );
    assert_eq!(
        primitives::re(&value),
        C::new(zero.from_i64(3), zero.clone())
    );
    assert_eq!(
        primitives::im(&value),
        C::new(zero.from_i64(-4), zero.clone())
    );
    assert_eq!(
        primitives::add(&value, &primitives::conj(&value)),
        C::new(zero.from_i64(6), zero.clone())
    );
    assert!(primitives::truth(&C::new(zero.clone(), zero.one())));
    assert_eq!(zero.from_i64(17).exact_i64(), Some(17));
    assert_eq!((zero.one() / zero.from_i64(2)).exact_i64(), None);
}

#[test]
fn native_primitives_exact_axes_all_types() {
    axes::<f64>();
    axes::<DoubleFloat>();
    axes::<Float>();
    assert_eq!(f64::FIXED_BITS, Some(53));
    assert_eq!(DoubleFloat::FIXED_BITS, Some(106));
    assert_eq!(Float::FIXED_BITS, None);
    assert_eq!(Float::DEFAULT_BITS, 128);
    assert_eq!(9_223_372_036_854_775_808.0_f64.exact_i64(), None);
    let large_integer = DoubleFloat::from_compensated_sum(2.0_f64.powi(60), 1.0);
    assert_eq!(large_integer.exact_i64(), Some((1_i64 << 60) + 1));
}

#[test]
fn native_primitives_binary64_hypot_exact_axes_and_nonfinite_values() {
    let values = [
        0.0,
        -0.0,
        f64::from_bits(1),
        -f64::from_bits(1),
        f64::MIN_POSITIVE,
        -f64::MIN_POSITIVE,
        1.0,
        -1.0,
        3.0,
        -4.0,
        f64::MAX,
        -f64::MAX,
        f64::INFINITY,
        f64::NEG_INFINITY,
        f64::NAN,
        -f64::NAN,
    ];
    for x in values {
        for y in values {
            let actual = NativeFloat::hypot(&x, &y);
            let expected = f64::hypot(x, y);
            if expected.is_nan() {
                assert!(actual.is_nan(), "hypot({x:?}, {y:?}) = {actual:?}");
            } else {
                assert_eq!(actual.to_bits(), expected.to_bits(), "hypot({x:?}, {y:?})");
            }
        }
    }
    for (x, y) in [(f64::NAN, f64::INFINITY), (f64::NEG_INFINITY, f64::NAN)] {
        assert_eq!(NativeFloat::hypot(&x, &y), f64::INFINITY);
    }
}

#[test]
fn native_primitives_disparate_component_precision() {
    const BITS: u32 = 512;
    let tiny = Float::parse("1e-150", Some(2)).unwrap();
    for real_sign in [-1, 1] {
        for imaginary_sign in [-1, 1] {
            let value = C::new(
                tiny.clone() * Float::with_val(2, real_sign),
                Float::with_val(BITS, 1.25 * f64::from(imaginary_sign)),
            );
            let exact = C::new(value.re.at_precision(1024), value.im.at_precision(1024));
            let square = primitives::powi(&value, 2, false);
            let expected_square = exact.clone() * &exact;
            // The real square remains well determined even though the tiny
            // input's relative precision is low. An exact unit at two bits
            // used to round -1.5625 to -1.5 in the final identity multiply.
            relative(&square.re, &expected_square.re, 480);
            assert!(square.re.prec() >= 480);
            let root = primitives::sqrt(&value);
            let expected_root = exact.sqrt();
            relative(&root.re, &expected_root.re, 480);
            relative(&root.im, &expected_root.im, 480);
            // Selecting a stronger prototype for exact constants must not
            // mutate inputs or pad a no-op result's component precision.
            assert_eq!(value.re.prec(), tiny.prec());
            let identity = primitives::powi(&value, 1, false);
            assert_eq!(identity.re.prec(), value.re.prec());
            assert_eq!(identity.im.prec(), value.im.prec());
            assert!(
                square.im.prec() <= tiny.prec() + 2,
                "square imaginary precision={}, input tiny precision={}",
                square.im.prec(),
                tiny.prec()
            );
        }
    }
    for sign in [-1, 1] {
        let axis = C::new(Float::with_val(2, -2), Float::with_val(BITS, 0.));
        let axis = if sign < 0 {
            C::new(axis.re, -axis.im)
        } else {
            axis
        };
        let result = primitives::log(&axis);
        let expected = Float::with_val(BITS, sign).pi() * Float::with_val(BITS, sign);
        relative(&result.im, &expected, 480);
        assert!(result.re.prec() < BITS);
        assert!(result.im.prec() >= BITS);

        let axis = C::new(Float::new(BITS), Float::with_val(2, sign));
        let result = primitives::log(&axis);
        relative(&result.im, &(expected / Float::with_val(BITS, 2)), 480);
        assert!(result.im.prec() >= BITS);
    }
}

#[test]
fn native_float_hypot_tracks_component_precision() {
    for bits in [512, 3456] {
        let tiny_text = if bits == 512 { "1e-170" } else { "1e-1100" };
        let tiny = Float::parse(tiny_text, Some(bits)).unwrap();
        let weak_tiny = Float::parse(tiny_text, Some(2)).unwrap();
        let weak_large = Float::with_val(2, 1.5);
        let strong_large = Float::with_val(bits, 1.5);
        let zero = Float::new(bits);
        for (large, small, maximum_bits, minimum_bits) in [
            (&weak_large, &zero, 2, 2),
            (&weak_large, &tiny, 2, 2),
            // Comparable components may gain a few metadata bits through
            // tracked addition/sqrt, but cannot inherit hundreds of bits
            // merely from the other input's representation.
            (&weak_large, &strong_large, 6, 2),
            (&strong_large, &weak_tiny, bits, bits - 16),
            (&strong_large, &tiny, bits, bits - 16),
            (&strong_large, &strong_large, bits, bits - 16),
        ] {
            for sign_x in [-1, 1] {
                for sign_y in [-1, 1] {
                    let x = large.clone() * Float::with_val(large.prec(), sign_x);
                    let y = small.clone() * Float::with_val(small.prec(), sign_y);
                    for (x, y) in [(&x, &y), (&y, &x)] {
                        let actual = NativeFloat::hypot(x, y);
                        assert!(
                            (minimum_bits..=maximum_bits).contains(&actual.prec()),
                            "hypot({x:?}, {y:?}) = {actual:?} at {} bits",
                            actual.prec()
                        );
                        let expected =
                            symbolica::domains::backend::float::MultiPrecisionFloat::with_val(
                                bits + 128,
                                x.as_raw().hypot_ref(y.as_raw()),
                            );
                        let delta =
                            symbolica::domains::backend::float::MultiPrecisionFloat::with_val(
                                bits + 128,
                                actual.as_raw() - &expected,
                            )
                            .abs();
                        let mut tolerance = expected;
                        tolerance >>= actual.prec().saturating_sub(4);
                        assert!(delta <= tolerance, "hypot({x:?}, {y:?}) = {actual:?}");
                    }
                }
            }
        }
        for x in [0., -0., f64::INFINITY, f64::NEG_INFINITY, f64::NAN] {
            for y in [0., -0., f64::INFINITY, f64::NEG_INFINITY, f64::NAN] {
                let x = Float::with_val(2, x);
                let y = Float::with_val(bits, y);
                let actual = NativeFloat::hypot(&x, &y);
                let expected = symbolica::domains::backend::float::MultiPrecisionFloat::with_val(
                    bits,
                    x.as_raw().hypot_ref(y.as_raw()),
                );
                if expected.is_nan() {
                    assert!(actual.as_raw().is_nan());
                } else {
                    assert_eq!(actual.as_raw(), &expected);
                    assert_eq!(actual.is_sign_negative(), expected.is_sign_negative());
                }
                assert_eq!(actual.prec(), bits);
            }
        }
    }
}

#[test]
fn native_primitives_binary64_subnormals_and_extremes() {
    for width in [f64::from_bits(1), 1e-300, -f64::from_bits(1), -1e-300] {
        assert!(primitives::truth(&C::new(0.0, width)));
        for real in [-1.0, 0.0, 1.0, width] {
            let value = C::new(real, width);
            let root = primitives::sqrt(&value);
            // The independent polar reference needs enough guard precision
            // to resolve a subnormal angle near the negative real axis.
            let exact = C::new(Float::with_val(2304, real), Float::with_val(2304, width)).sqrt();
            for (actual, expected) in [(root.re, exact.re.to_f64()), (root.im, exact.im.to_f64())] {
                assert!(
                    (actual - expected).abs()
                        <= 8.0 * f64::EPSILON * expected.abs() + f64::from_bits(2),
                    "{value:?}: {root:?} vs {exact:?}"
                );
            }
            assert_eq!(root.im.is_sign_negative(), width.is_sign_negative());
        }
    }
    for value in [
        C::new(f64::MAX, f64::MAX),
        C::new(-f64::MAX, f64::MAX),
        C::new(1e-300, -2e-300),
    ] {
        let root = primitives::sqrt(&value);
        let logarithm = primitives::log(&value);
        let reciprocal = primitives::powi(&value, -1, false);
        let exact = C::new(
            Float::with_val(256, value.re),
            Float::with_val(256, value.im),
        );
        for (actual, expected) in [
            (root, exact.sqrt()),
            (logarithm, exact.log()),
            (reciprocal, exact.inv()),
        ] {
            for (actual, expected) in [
                (actual.re, expected.re.to_f64()),
                (actual.im, expected.im.to_f64()),
            ] {
                assert!(actual.is_finite());
                assert!(
                    (actual - expected).abs()
                        <= 16.0 * f64::EPSILON * expected.abs() + f64::from_bits(2)
                );
            }
        }
    }
    let value = C::new(0.5717239802032432, -0.06483875997714851);
    assert_eq!(primitives::mul(&value, &primitives::conj(&value)).im, 0.0);
}

#[test]
fn native_primitives_double_float_retains_low_components() {
    let value = C::new(
        DoubleFloat::from_compensated_sum(0.5, 1e-20),
        DoubleFloat::from_compensated_sum(0.25, -2e-20),
    );
    let reference = C::new(dd_float(value.re, 256), dd_float(value.im, 256));
    for (actual, expected) in [
        (primitives::sqrt(&value), reference.sqrt()),
        (primitives::log(&value), reference.log()),
        (primitives::powi(&value, -1, false), reference.inv()),
        (primitives::polylog2(&value), Float::dilog(&reference)),
    ] {
        relative(&dd_float(actual.re, 256), &expected.re, 106);
        relative(&dd_float(actual.im, 256), &expected.im, 106);
    }
    let rounded = C::new(DoubleFloat::from(0.5), DoubleFloat::from(0.25));
    assert_ne!(primitives::polylog2(&value), primitives::polylog2(&rounded));
    for width in [-1e-300, 1e-300] {
        let lower = C::new(DoubleFloat::from(-1.0), DoubleFloat::from(width));
        let root = primitives::sqrt(&lower);
        assert!(primitives::truth(&lower));
        assert_eq!(root.im.negative_sign(), width.is_sign_negative());
        assert_eq!(root.re, DoubleFloat::from(width.abs() / 2.0));
    }
}

#[test]
fn native_primitives_double_float_products_keep_exact_symmetries() {
    // The compensated DoubleFloat product rounds its cross terms in operand
    // order. The native complex product forms them symmetrically, so
    // `z conj(z)` is exactly real (as for binary64 above) and products
    // commute. A residue of 1e-35 in `1 - z conj(z)` made a four-mass D0
    // invert a dilogarithm argument of modulus 0.09
    // (tests/native_double_float_d0.rs).
    let ratio = |n: i64, d: i64| DoubleFloat::from(n as f64) / DoubleFloat::from(d as f64);
    let mut residues = 0;
    for n in 1..=40 {
        let z = C::new(ratio(n, 123), ratio(41 - n, 287));
        let w = C::new(ratio(n + 3, 11), -ratio(n, 13));
        let norm = primitives::mul(&z, &primitives::conj(&z));
        assert!(norm.im.is_fully_zero(), "{z:?} conj: {norm:?}");
        let negated = primitives::mul(&C::new(-z.re, -z.im), &primitives::conj(&z));
        assert!(
            negated.im.is_fully_zero(),
            "{z:?} negated conj: {negated:?}"
        );
        assert_eq!(primitives::mul(&z, &w), primitives::mul(&w, &z));
        if !(z.re * -z.im + z.im * z.re).is_fully_zero() {
            residues += 1;
        }
    }
    // The raw products leave such residues for some of these values.
    assert!(residues > 0);
}

#[test]
fn native_primitives_general_powers_and_quadrants() {
    for re in [-2.0, -0.25, 0.25, 2.0] {
        for im in [-0.75, 0.75] {
            let value = C::new(re, im);
            let reference = C::new(Float::with_val(256, re), Float::with_val(256, im));
            for exponent in [C::new(0.25, 0.0), C::new(-0.25, 0.125), C::new(3.0, 0.0)] {
                let actual = primitives::powf(&value, &exponent, false);
                let expected = reference.powf(&C::new(
                    Float::with_val(256, exponent.re),
                    Float::with_val(256, exponent.im),
                ));
                for (a, b) in [
                    (actual.re, expected.re.to_f64()),
                    (actual.im, expected.im.to_f64()),
                ] {
                    assert!((a - b).abs() <= 32.0 * f64::EPSILON * b.abs());
                }
                let value = C::new(DoubleFloat::from(re), DoubleFloat::from(im));
                let exponent = C::new(
                    DoubleFloat::from(exponent.re),
                    DoubleFloat::from(exponent.im),
                );
                let actual = primitives::powf(&value, &exponent, false);
                relative(&dd_float(actual.re, 256), &expected.re, 106);
                relative(&dd_float(actual.im, 256), &expected.im, 106);
            }
            for exponent in [-3, -1, 0, 1, 2, 3] {
                let actual = primitives::powi(&value, exponent, false);
                let expected =
                    reference.powf(&C::new(Float::with_val(256, exponent), Float::new(256)));
                for (a, b) in [
                    (actual.re, expected.re.to_f64()),
                    (actual.im, expected.im.to_f64()),
                ] {
                    assert!((a - b).abs() <= 32.0 * f64::EPSILON * (1.0 + b.abs()));
                }
            }
        }
    }
    for zero in [0.0, -0.0] {
        let actual = primitives::powf(&C::new(-4.0, zero), &C::new(0.25, 0.0), true);
        assert!((actual.re - 1.0).abs() < 4.0 * f64::EPSILON);
        assert!((actual.im.abs() - 1.0).abs() < 4.0 * f64::EPSILON);
        assert_eq!(actual.im.is_sign_negative(), zero.is_sign_negative());
        let dilog = primitives::polylog2(&C::new(2.0, zero));
        // The registered evaluator rounds once to binary64; the analytic
        // reference here rounds pi, log(2) and their product separately.
        let expected = -std::f64::consts::PI * 2.0_f64.ln();
        assert!((dilog.im - expected).abs() <= 4.0 * f64::EPSILON * expected.abs());
    }
}

#[test]
fn native_binary64_dilog_known_real_and_complex_points() {
    let pi = std::f64::consts::PI;
    let log_two = 2.0_f64.ln();
    let catalan = 0.915_965_594_177_219;
    for (input, expected) in [
        (C::new(0.0, 0.0), C::new(0.0, 0.0)),
        (C::new(1.0, 0.0), C::new(pi * pi / 6.0, 0.0)),
        (C::new(-1.0, 0.0), C::new(-pi * pi / 12.0, 0.0)),
        (
            C::new(0.5, 0.0),
            C::new(pi * pi / 12.0 - log_two * log_two / 2.0, 0.0),
        ),
        (C::new(0.0, 1.0), C::new(-pi * pi / 48.0, catalan)),
        (C::new(0.0, -1.0), C::new(-pi * pi / 48.0, -catalan)),
        (
            C::new(1.0, 1.0),
            C::new(pi * pi / 16.0, catalan + pi * log_two / 4.0),
        ),
    ] {
        let actual = primitives::polylog2(&input);
        for (value, reference) in [(actual.re, expected.re), (actual.im, expected.im)] {
            assert!(
                (value - reference).abs() <= 8.0 * f64::EPSILON * reference.abs().max(1.0),
                "Li2({input:?}) = {actual:?}, expected {expected:?}"
            );
        }
        if expected.im == 0.0 {
            assert_eq!(actual.im, 0.0);
        }
    }
}

#[test]
fn native_binary64_dilog_cut_sides_and_signed_zero() {
    let pi = std::f64::consts::PI;
    let expected_real = pi * pi / 4.0;
    let imaginary_magnitude = pi * 2.0_f64.ln();
    // Exactly real inputs take the lower lip for BOTH signs of zero.
    // Nonzero imaginary parts instead approach the corresponding cut side.
    for imaginary in [0.0, -0.0, 1e-20, -1e-20, 1e-300, -1e-300] {
        let actual = primitives::polylog2(&C::new(2.0, imaginary));
        let expected_imaginary = if imaginary > 0.0 {
            imaginary_magnitude
        } else {
            -imaginary_magnitude
        };
        assert!((actual.re - expected_real).abs() <= 8.0 * f64::EPSILON * expected_real);
        assert!(
            (actual.im - expected_imaginary).abs()
                <= 8.0 * f64::EPSILON * imaginary_magnitude,
            "Li2(2 + {imaginary:?} i) = {actual:?}"
        );
        assert_eq!(actual.im.is_sign_negative(), imaginary <= 0.0);
    }
}

#[test]
fn native_scalar_half_powers_preserve_unbalanced_precision() {
    for (bits, tiny, tolerance) in [(512, "1e-170", "1e-140"), (3456, "1e-1100", "1e-1000")] {
        let n = |text| Float::parse(text, Some(bits)).unwrap();
        let tiny = Float::parse(tiny, Some(2)).unwrap();
        for real_negative in [false, true] {
            for imaginary_negative in [false, true] {
                for tiny_real in [false, true] {
                    let re = if tiny_real { tiny.clone() } else { n("2") };
                    let im = if tiny_real { n("1.25") } else { tiny.clone() };
                    let value = C::new(
                        if real_negative { -re } else { re },
                        if imaginary_negative { -im } else { im },
                    );
                    let refined = C::new(
                        value.re.at_precision(bits + 128),
                        value.im.at_precision(bits + 128),
                    );
                    let root = refined.sqrt();
                    for exponent in ["-1.5", "-0.5", "0.5", "1.5"] {
                        let actual = primitives::powf(&value, &C::new(n(exponent), n("0")), false);
                        let expected = if exponent.ends_with("1.5") {
                            refined.clone() * root.clone()
                        } else {
                            root.clone()
                        };
                        let expected = if exponent.starts_with('-') {
                            expected.inv()
                        } else {
                            expected
                        };
                        assert!(
                            (actual.re.clone() - expected.re).norm() < n(tolerance),
                            "exponent={exponent}"
                        );
                        assert!(
                            (actual.im.clone() - expected.im).norm() < n(tolerance),
                            "exponent={exponent}"
                        );
                    }
                }
            }
        }
    }
}

#[test]
fn native_inverse_three_halves_avoids_overflowing_intermediates() {
    for re in [1e210, -1e210, 1e200, -1e200] {
        for im in [0.0, -0.0, re * 0.5, -re * 0.5] {
            let value = C::new(re, im);
            let actual = primitives::powf(&value, &C::new(-1.5, 0.0), false);
            let refined = C::new(Float::with_val(256, re), Float::with_val(256, im));
            let expected = refined.powf(&C::new(Float::with_val(256, -1.5), Float::new(256)));
            for (a, e) in [
                (actual.re, expected.re.to_f64()),
                (actual.im, expected.im.to_f64()),
            ] {
                assert!(a.is_finite());
                assert!(
                    (a - e).abs() <= 64.0 * f64::EPSILON * e.abs() + f64::from_bits(16),
                    "a={a:e}, expected={e:e}"
                );
                if e.abs() > f64::from_bits(16) {
                    assert_ne!(a, 0.0);
                }
            }
        }
    }
}

#[test]
fn native_integer_powers_preserve_real_axis_signed_lips() {
    fn check<T: NativeFloat>() {
        let n = T::zero_at(512);
        for real in [-2, 2] {
            for imaginary_negative in [false, true] {
                let zero = n.zero();
                let z = C::new(
                    n.from_i64(real),
                    if imaginary_negative { -zero } else { zero },
                );
                for exponent in [-5, -4, -3, -2, -1, 1, 2, 3, 4, 5] {
                    let actual = primitives::powi(&z, exponent, false);
                    assert!(actual.im.is_fully_zero());
                    let expected_negative =
                        imaginary_negative ^ (exponent < 0) ^ (real < 0 && exponent % 2 == 0);
                    assert_eq!(
                        actual.im.negative_sign(),
                        expected_negative,
                        "type={}, base={real}, exponent={exponent}, imaginary_negative={imaginary_negative}",
                        std::any::type_name::<T>()
                    );
                }
            }
        }
    }
    check::<f64>();
    // DoubleFloat's compensated arithmetic has independent signed-zero
    // normalization. These non-inverse integer powers are not present in the
    // current scalar DAG (its integer-power instructions are all -1).
    check::<Float>();
}

#[test]
fn native_primitives_thousand_decimal_digits() {
    const BITS: u32 = 3456;
    let value = C::new(Float::with_val(BITS, 0.25), Float::with_val(BITS, 0.5));
    let root = primitives::sqrt(&value);
    let squared = primitives::mul(&root, &root);
    relative(&squared.re, &value.re, BITS);
    relative(&squared.im, &value.im, BITS);
    let actual = primitives::polylog2(&value);
    let refined = C::new(
        value.re.at_precision(BITS + 384),
        value.im.at_precision(BITS + 384),
    );
    let expected = Float::dilog(&refined);
    relative(&actual.re, &expected.re, BITS);
    relative(&actual.im, &expected.im, BITS);
    let width = Float::parse("-1e-1000", Some(BITS)).unwrap();
    let lower = C::new(Float::with_val(BITS, -1), width.clone());
    let root = primitives::sqrt(&lower);
    assert!(root.im.negative_sign());
    // Re sqrt(-1+iy)=|y|/2+O(y^3), checked against the tiny component itself.
    relative(&root.re, &(width.norm() / Float::with_val(BITS, 2)), BITS);
    assert!(primitives::truth(&C::new(Float::new(BITS), width)));
}

/// Keep exceptional inputs visible while retaining the finite-axis optimization.
fn reciprocal_axes<T: NativeFloat>(from_f64: impl Fn(f64) -> T) {
    let zero = from_f64(0.0);
    for value in [1.0, -1.0, 1e110, -1e110, 1e-110, -1e-110] {
        let value = from_f64(value);
        for axis in [
            C::new(value.clone(), zero.clone()),
            C::new(zero.clone(), value.clone()),
        ] {
            let actual = primitives::powi(&axis, -1, false);
            assert!(actual.is_finite());
            if axis.im.is_fully_zero() {
                assert_eq!(actual.re, value.inv());
                assert!(actual.im.is_fully_zero());
            } else {
                assert!(actual.re.is_fully_zero());
                assert_eq!(actual.im, -value.inv());
            }
        }
    }
    for real in [0.0, -0.0] {
        for imaginary in [0.0, -0.0] {
            assert!(
                !primitives::powi(&C::new(from_f64(real), from_f64(imaginary)), -1, false)
                    .is_finite()
            );
        }
    }
    for real in [0.0, -0.0, 1.0, f64::INFINITY, f64::NEG_INFINITY, f64::NAN] {
        for imaginary in [0.0, -0.0, 1.0, f64::INFINITY, f64::NEG_INFINITY, f64::NAN] {
            if real.is_finite() && imaginary.is_finite() {
                continue;
            }
            let input = C::new(from_f64(real), from_f64(imaginary));
            let actual = primitives::powi(&input, -1, false);
            let expected = input.inv();
            assert!(
                !actual.is_finite(),
                "nonfinite input became finite: {input:?} -> {actual:?}"
            );
            for (actual, expected) in [(&actual.re, &expected.re), (&actual.im, &expected.im)] {
                if expected.real_classify() == Some(std::num::FpCategory::Nan) {
                    assert_eq!(actual.real_classify(), Some(std::num::FpCategory::Nan));
                } else {
                    assert_eq!(actual, expected);
                    assert_eq!(actual.negative_sign(), expected.negative_sign());
                }
                assert_eq!(actual.get_precision(), expected.get_precision());
            }
        }
    }
}

#[test]
fn native_reciprocal_nonfinite_axes_delegate_to_symbolica() {
    reciprocal_axes::<f64>(|x| x);
    reciprocal_axes::<DoubleFloat>(DoubleFloat::from);
    reciprocal_axes::<Float>(|x| Float::with_val(256, x));
}

#[test]
fn native_binary64_reciprocal_and_norm_fast_paths_stay_accurate() {
    // Off-axis binary64 values in [2^-500, 2^500] use one division for the
    // reciprocal and sqrt(x²+y²) for the norm; other values keep the scaled
    // algorithms. Both must stay within a few ulps of a 256-bit reference.
    const BITS: u32 = 256;
    let reference_inverse = |z: C<f64>| {
        let (x, y) = (Float::with_val(BITS, z.re), Float::with_val(BITS, z.im));
        let norm = x.clone() * &x + y.clone() * &y;
        (x / &norm, -y / norm)
    };
    let mut state = 0x9e37_79b9_7f4a_7c15_u64;
    let mut next = || {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        (state >> 11) as f64 / (1u64 << 53) as f64
    };
    let mut points = Vec::new();
    for _ in 0..20_000 {
        let magnitude = 2f64.powf(1200. * next() - 600.);
        let ratio = 2f64.powf(120. * next() - 60.);
        let signs = [(1., 1.), (-1., 1.), (1., -1.), (-1., -1.)][(4. * next()) as usize % 4];
        points.push(C::new(signs.0 * magnitude, signs.1 * magnitude * ratio));
    }
    points.extend([
        C::new(1e300, 1e300),
        C::new(-1e-300, 1e-300),
        C::new(3e150, 1.0),
        C::new(3.1e-151, -2e-160),
    ]);
    for z in points {
        let inverse = primitives::inverse(&z);
        let (re, im) = reference_inverse(z);
        let scale = re.to_f64().abs().max(im.to_f64().abs());
        for (value, expected) in [(inverse.re, re), (inverse.im, im)] {
            let error = (Float::with_val(BITS, value) - expected).norm().to_f64();
            assert!(error <= 4.0 * f64::EPSILON * scale, "1/{z:?} = {inverse:?}");
        }
        let norm = NativeFloat::hypot(&z.re, &z.im);
        let expected =
            (Float::with_val(BITS, z.re).pow(2) + Float::with_val(BITS, z.im).pow(2)).sqrt();
        let error = (Float::with_val(BITS, norm) - &expected).norm().to_f64();
        assert!(
            error <= 2.0 * f64::EPSILON * expected.to_f64(),
            "|{z:?}| = {norm:e}"
        );
    }
}
