//! Warm direct-Rust throughput over mixed master rows, matching
//! tests/support/rows_bench_oracle.f90 (see tests/generate_bench_rows.py).
//! Usage: rows_survey ROWS [calls=20000] [repetitions=7] [batches=1,1024] [--values]
//! Kinds 1-4 use `NativeEvaluator<f64>` (A0, B0, C0, D0); 5-7 use
//! `NativeTensorEvaluator<f64>` (An, B11, B111). Workloads: SAME repeats a row,
//! HETERO cycles one kind's rows, ALL cycles every row in file order. Batch 1
//! evaluates each row of the stream in order; batch B > 1 evaluates each kind's
//! rows of the same stream with `evaluate_batch` in chunks of B rows. Input
//! grouping, allocation and validation of outputs are outside timing; per-row
//! input validation inside the public API is included. A timed block has at
//! least `calls` evaluations, raised after an untimed calibration block so that
//! it lasts at least 25 ms. Emits TSV.
use oneloop::{NativeEvaluator, NativeTensorEvaluator, ScalarIntegral, TensorIntegral};
use std::{hint::black_box, time::Instant};
use symbolica::domains::float::Complex;

type C = Complex<f64>;

struct Row {
    name: String,
    kind: usize,
    input: Vec<C>,
}

enum Evaluator {
    Scalar(NativeEvaluator<f64>),
    Tensor(NativeTensorEvaluator<f64>),
}

impl Evaluator {
    fn new(kind: usize) -> Self {
        match kind {
            1 => Self::Scalar(NativeEvaluator::new(ScalarIntegral::A0).unwrap()),
            2 => Self::Scalar(NativeEvaluator::new(ScalarIntegral::B0).unwrap()),
            3 => Self::Scalar(NativeEvaluator::new(ScalarIntegral::C0).unwrap()),
            4 => Self::Scalar(NativeEvaluator::new(ScalarIntegral::D0).unwrap()),
            5 => Self::Tensor(NativeTensorEvaluator::new(TensorIntegral::An).unwrap()),
            6 => Self::Tensor(NativeTensorEvaluator::new(TensorIntegral::B11).unwrap()),
            7 => Self::Tensor(NativeTensorEvaluator::new(TensorIntegral::B111).unwrap()),
            _ => panic!("unsupported row kind {kind}"),
        }
    }

    fn outputs(&self) -> usize {
        match self {
            Self::Scalar(_) => 3,
            Self::Tensor(e) => e.family().outputs(),
        }
    }

    fn arity(&self) -> usize {
        match self {
            Self::Scalar(e) => e.family().arity(),
            Self::Tensor(e) => e.family().arity(),
        }
    }

    #[inline]
    fn evaluate(&mut self, input: &[C], output: &mut [C]) {
        match self {
            Self::Scalar(e) => e.evaluate(input, output).unwrap(),
            Self::Tensor(e) => e.evaluate(input, output).unwrap(),
        }
    }

    #[inline]
    fn evaluate_batch(&mut self, input: &[C], output: &mut [C], rows: usize) {
        match self {
            Self::Scalar(e) => e.evaluate_batch(input, output, rows).unwrap(),
            Self::Tensor(e) => e.evaluate_batch(input, output, rows).unwrap(),
        }
    }
}

fn parse(path: &str) -> Vec<Row> {
    let mut rows = Vec::new();
    let mut name = String::new();
    for line in std::fs::read_to_string(path).unwrap().lines() {
        if let Some(rest) = line.strip_prefix("# name ") {
            name = rest.trim().into();
            continue;
        }
        if line.starts_with('#') || line.trim().is_empty() {
            continue;
        }
        let values: Vec<f64> = line
            .split_whitespace()
            .map(|v| v.parse().unwrap())
            .collect();
        let kind = values[0] as usize;
        let (momenta, masses) = match kind {
            1 | 5 => (0, 1),
            2 | 6 | 7 => (1, 2),
            3 => (3, 3),
            4 => (6, 4),
            _ => panic!("unsupported row kind {kind}"),
        };
        assert_eq!(values.len(), 2 + momenta + 2 * masses, "{line}");
        let mut input: Vec<C> = values[2..2 + momenta]
            .iter()
            .map(|&p| C::new(p, 0.))
            .collect();
        input.extend(
            values[2 + momenta..]
                .chunks_exact(2)
                .map(|m| C::new(m[0], m[1])),
        );
        input.push(C::new(values[1], 0.));
        rows.push(Row {
            name: if name.is_empty() {
                format!("row_{}", rows.len() + 1)
            } else {
                std::mem::take(&mut name)
            },
            kind,
            input,
        });
    }
    rows
}

/// One timed stream: either row by row, or per-kind batches of `batch` rows.
struct Stream {
    order: Vec<usize>,
    groups: Vec<(usize, Vec<C>, Vec<C>, usize)>,
}

impl Stream {
    fn new(
        rows: &[Row],
        selection: &[usize],
        calls: usize,
        evaluators: &[Option<Evaluator>],
    ) -> Self {
        let order: Vec<usize> = selection.iter().copied().cycle().take(calls).collect();
        let mut groups = Vec::new();
        for kind in 1..=7 {
            let Some(evaluator) = &evaluators[kind] else {
                continue;
            };
            let members: Vec<_> = order
                .iter()
                .copied()
                .filter(|&r| rows[r].kind == kind)
                .collect();
            if members.is_empty() {
                continue;
            }
            let input: Vec<C> = members
                .iter()
                .flat_map(|&r| rows[r].input.iter().copied())
                .collect();
            let output = vec![C::new(f64::NAN, f64::NAN); members.len() * evaluator.outputs()];
            groups.push((kind, input, output, members.len()));
        }
        Self { order, groups }
    }

    fn run(
        &mut self,
        rows: &[Row],
        evaluators: &mut [Option<Evaluator>],
        batch: usize,
        scratch: &mut [C],
    ) -> C {
        let mut checksum = C::new(0., 0.);
        if batch == 1 {
            for &r in &self.order {
                let evaluator = evaluators[rows[r].kind].as_mut().unwrap();
                let width = evaluator.outputs();
                evaluator.evaluate(black_box(&rows[r].input), &mut scratch[..width]);
                checksum += scratch[0];
            }
        } else {
            for (kind, input, output, count) in &mut self.groups {
                let evaluator = evaluators[*kind].as_mut().unwrap();
                let (arity, width) = (evaluator.arity(), evaluator.outputs());
                for (chunk_in, chunk_out) in input
                    .chunks(arity * batch)
                    .zip(output.chunks_mut(width * batch))
                {
                    evaluator.evaluate_batch(
                        black_box(chunk_in),
                        chunk_out,
                        chunk_in.len() / arity,
                    );
                }
                for row in output.chunks_exact(width).take(*count) {
                    checksum += row[0];
                }
            }
        }
        checksum
    }
}

fn main() {
    let mut args: Vec<String> = std::env::args().skip(1).collect();
    let values = args.iter().any(|a| a == "--values");
    args.retain(|a| a != "--values");
    let path = args.first().expect("rows file").clone();
    let calls: usize = args.get(1).map_or(20000, |s| s.parse().unwrap());
    let repetitions: usize = args.get(2).map_or(7, |s| s.parse().unwrap());
    let batches: Vec<usize> = args
        .get(3)
        .map_or("1,1024".into(), |s| s.clone())
        .split(',')
        .map(|s| s.parse().unwrap())
        .collect();
    std::thread::Builder::new()
        .stack_size(128 * 1024 * 1024)
        .spawn(move || run(&path, calls, repetitions, &batches, values))
        .unwrap()
        .join()
        .unwrap();
}

fn run(path: &str, calls: usize, repetitions: usize, batches: &[usize], values: bool) {
    let rows = parse(path);
    let mut evaluators: Vec<Option<Evaluator>> = (0..=7).map(|_| None).collect();
    for row in &rows {
        if evaluators[row.kind].is_none() {
            evaluators[row.kind] = Some(Evaluator::new(row.kind));
        }
    }
    let mut scratch = vec![C::new(0., 0.); 12];
    if values {
        for (index, row) in rows.iter().enumerate() {
            let evaluator = evaluators[row.kind].as_mut().unwrap();
            let width = evaluator.outputs();
            evaluator.evaluate(&row.input, &mut scratch[..width]);
            print!("VALUE\t{}\t{}", index + 1, row.kind);
            for z in &scratch[..width] {
                print!("\t{:.17e}\t{:.17e}", z.re, z.im);
            }
            println!();
        }
    }
    println!("mode\tkind\tname\tbatch\tcalls\trepetition\tns_per_call\tchecksum_re\tchecksum_im");
    let mut workloads: Vec<(String, usize, String, Vec<usize>)> = Vec::new();
    for (index, row) in rows.iter().enumerate() {
        workloads.push(("SAME".into(), row.kind, row.name.clone(), vec![index]));
    }
    for kind in 1..=7 {
        let selection: Vec<_> = (0..rows.len()).filter(|&r| rows[r].kind == kind).collect();
        if !selection.is_empty() {
            workloads.push(("HETERO".into(), kind, "all".into(), selection));
        }
    }
    workloads.push(("ALL".into(), 0, "all".into(), (0..rows.len()).collect()));
    for &batch in batches {
        for (mode, kind, name, selection) in &workloads {
            let mut warm = Stream::new(&rows, selection, 256.min(calls), &evaluators);
            warm.run(&rows, &mut evaluators, batch, &mut scratch);
            let mut stream = Stream::new(&rows, selection, calls, &evaluators);
            let started = Instant::now();
            stream.run(&rows, &mut evaluators, batch, &mut scratch);
            let seconds = started.elapsed().as_secs_f64();
            let block = if seconds < 0.025 {
                (calls as f64 * 0.025 / seconds.max(1e-6)) as usize + 1
            } else {
                calls
            };
            if block != calls {
                stream = Stream::new(&rows, selection, block, &evaluators);
            }
            for repetition in 1..=repetitions {
                let started = Instant::now();
                let checksum = stream.run(&rows, &mut evaluators, batch, &mut scratch);
                let elapsed = started.elapsed().as_nanos() as f64 / block as f64;
                println!(
                    "{mode}\t{kind}\t{name}\t{batch}\t{block}\t{repetition}\t{elapsed:.6}\t{:.17e}\t{:.17e}",
                    checksum.re, checksum.im
                );
            }
        }
    }
}
