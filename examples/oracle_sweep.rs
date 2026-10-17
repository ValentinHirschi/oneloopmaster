//! Evaluate sweep rows at binary64 and high precision.
//!
//! Run: cargo run --release --features generated-evaluators --example oracle_sweep
//!        [-- --bits 256 --backends --expression] < rows.txt
//! Rows use the scalar fixture order without references:
//! `kind mu_squared momenta... (squared_mass_re squared_mass_im)...`, with
//! kinds 1,2,-2,3,4 for A0,B0,dB0,C0,D0. Each row prints `f64` and `fN` lines
//! with `[finite, 1/eps, 1/eps^2]` real/imaginary pairs; `--backends` adds the
//! binary64 `expr` and `jit` lines of the Expression and SymJIT backends.
//! The `f64`/`fN` lines use the native backend when it is compiled in, and
//! otherwise (or with `--expression`) the exact-expression interpreter, so
//! expression changes can be swept before regenerating native code.
//! `--generated` evaluates binary64 with the generated code in every sector
//! instead of the hand-written kernels.
//! tests/oracle_sweep.py drives this together with tests/support/sweep_oracle.f90.
use oneloop::{EvaluationBackend, PrecisionEvaluator, ScalarIntegral, evaluate_with_backend};
use std::io::{BufRead, Write};
use symbolica::domains::float::{Complex, Float, RealLike};

type C = Complex<f64>;

/// Binary64 and high-precision evaluation of one family.
struct Lanes {
    family: ScalarIntegral,
    #[cfg(feature = "generated-evaluators")]
    native: Option<oneloop::NativeEvaluator<f64>>,
    precise: PrecisionEvaluator,
}

impl Lanes {
    fn new(family: ScalarIntegral, bits: u32, expression: bool, generated: bool) -> Self {
        #[cfg(not(feature = "generated-evaluators"))]
        let _ = generated;
        let backend = if expression || !EvaluationBackend::Native.is_available() {
            EvaluationBackend::Expression
        } else {
            EvaluationBackend::Native
        };
        Self {
            family,
            #[cfg(feature = "generated-evaluators")]
            native: (backend == EvaluationBackend::Native).then(|| {
                let mut native = oneloop::NativeEvaluator::<f64>::new(family).unwrap();
                native.set_handwritten(!generated);
                native
            }),
            precise: PrecisionEvaluator::with_binary_precision_and_backend(family, bits, backend)
                .unwrap(),
        }
    }

    fn binary64(&mut self, input: &[C], output: &mut [C]) -> Result<(), String> {
        #[cfg(feature = "generated-evaluators")]
        if let Some(native) = &mut self.native {
            return native.evaluate(input, output);
        }
        evaluate_with_backend(self.family, input, output, EvaluationBackend::Expression)
    }
}

fn run() {
    let mut bits = 256u32;
    let mut backends = false;
    let mut expression = false;
    let mut generated = false;
    let mut arguments = std::env::args().skip(1);
    while let Some(argument) = arguments.next() {
        match argument.as_str() {
            "--bits" => {
                bits = arguments
                    .next()
                    .and_then(|b| b.parse().ok())
                    .expect("--bits N")
            }
            "--backends" => backends = true,
            "--expression" => expression = true,
            "--generated" => generated = true,
            other => panic!("unknown argument {other}"),
        }
    }
    let families = [
        ScalarIntegral::A0,
        ScalarIntegral::B0,
        ScalarIntegral::DB0,
        ScalarIntegral::C0,
        ScalarIntegral::D0,
    ];
    // Prepare families lazily: expression interpreters are costly to build.
    let mut lanes: [Option<Lanes>; 5] = Default::default();
    let stdout = std::io::stdout();
    let mut out = stdout.lock();
    let format = |values: &[(f64, f64)]| {
        values
            .iter()
            .map(|(re, im)| format!("{re:.17e} {im:.17e}"))
            .collect::<Vec<_>>()
            .join(" ")
    };
    for line in std::io::stdin().lock().lines() {
        let line = line.unwrap();
        let values = line
            .split_whitespace()
            .map(|v| v.parse::<f64>().unwrap())
            .collect::<Vec<_>>();
        if values.is_empty() {
            continue;
        }
        let (index, momenta, masses) = match values[0] as i64 {
            1 => (0, 0, 1),
            2 => (1, 1, 2),
            -2 => (2, 1, 2),
            3 => (3, 3, 3),
            4 => (4, 6, 4),
            kind => panic!("unsupported kind {kind}"),
        };
        assert_eq!(
            values.len(),
            2 + momenta + 2 * masses,
            "malformed row {line}"
        );
        let mut args = values[2..2 + momenta]
            .iter()
            .map(|&v| (v, 0.))
            .collect::<Vec<_>>();
        args.extend(values[2 + momenta..].chunks_exact(2).map(|v| (v[0], v[1])));
        args.push((values[1], 0.));

        let lane = lanes[index]
            .get_or_insert_with(|| Lanes::new(families[index], bits, expression, generated));
        let input = args
            .iter()
            .map(|&(re, im)| Complex::new(re, im))
            .collect::<Vec<_>>();
        let mut output = [Complex::new(f64::NAN, f64::NAN); 3];
        let result = lane.binary64(&input, &mut output);
        let row = output.map(|z| (z.re, z.im));
        match result {
            Ok(()) => writeln!(out, "f64 {}", format(&row)).unwrap(),
            Err(error) => writeln!(out, "f64 error {error}").unwrap(),
        }

        let input = args
            .iter()
            .map(|&(re, im)| Complex::new(Float::with_val(bits, re), Float::with_val(bits, im)))
            .collect::<Vec<_>>();
        let mut output =
            std::array::from_fn::<_, 3, _>(|_| Complex::new(Float::new(bits), Float::new(bits)));
        let result = lane.precise.evaluate(&input, &mut output);
        let row = output.map(|z| (z.re.to_f64(), z.im.to_f64()));
        match result {
            Ok(()) => writeln!(out, "f{bits} {}", format(&row)).unwrap(),
            Err(error) => writeln!(out, "f{bits} error {error}").unwrap(),
        }

        if backends {
            let input = args
                .iter()
                .map(|&(re, im)| Complex::new(re, im))
                .collect::<Vec<_>>();
            for (name, backend) in [
                ("expr", EvaluationBackend::Expression),
                ("jit", EvaluationBackend::SymJit),
            ] {
                let mut output = [Complex::new(f64::NAN, f64::NAN); 3];
                match evaluate_with_backend(families[index], &input, &mut output, backend) {
                    Ok(()) => {
                        writeln!(out, "{name} {}", format(&output.map(|z| (z.re, z.im)))).unwrap()
                    }
                    Err(error) => writeln!(out, "{name} error {error}").unwrap(),
                }
            }
        }
        out.flush().unwrap();
    }
}

fn main() {
    std::thread::Builder::new()
        .stack_size(128 * 1024 * 1024)
        .spawn(run)
        .unwrap()
        .join()
        .unwrap();
}
