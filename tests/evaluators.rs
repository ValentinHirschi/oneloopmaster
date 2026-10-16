//! Portable fresh-load and heterogeneous batch regressions against original data.
#[path = "support/fixtures.rs"]
mod fixtures;
use oneloop::{ScalarEvaluator, ScalarIntegral};
use symbolica::domains::float::Complex;

#[test]
fn eager_initialization_compiles_every_shared_master_body() {
    std::thread::Builder::new()
        .stack_size(128 * 1024 * 1024)
        .spawn(|| {
            oneloop::initialize().unwrap();
            assert!(oneloop::is_initialized());
        })
        .unwrap()
        .join()
        .unwrap();
}

#[test]
fn portable_scalar_caches_and_mixed_batches_match_original_fixtures() {
    std::thread::Builder::new()
        .stack_size(128 * 1024 * 1024)
        .spawn(|| {
            let fixtures = fixtures::parse(include_str!("data/parity.txt"));
            let mut failures = Vec::new();
            for family in [
                ScalarIntegral::A0,
                ScalarIntegral::B0,
                ScalarIntegral::DB0,
                ScalarIntegral::C0,
                ScalarIntegral::D0,
            ] {
                #[cfg(feature = "prebuilt")]
                let evaluator = ScalarEvaluator::prebuilt(family).unwrap();
                #[cfg(not(feature = "prebuilt"))]
                let evaluator = ScalarEvaluator::rebuild(family).unwrap();
                let bytes = evaluator.to_bytes().unwrap();
                drop(evaluator);
                let mut restored = ScalarEvaluator::from_bytes(family, &bytes).unwrap();
                assert!(ScalarEvaluator::from_bytes(family, &bytes[..1]).is_err());
                // Corrupt each size independently without touching valid executable
                // IR. Reject it before any undersized slice reaches generated code.
                // Skip the family discriminator, then the NUL-terminated version.
                let header_end = bytes[1..].iter().position(|&b| b == 0).unwrap() + 2;
                for offset in [header_end, header_end + 8] {
                    let mut corrupt = bytes.clone();
                    corrupt[offset..offset + 8].copy_from_slice(&0u64.to_le_bytes());
                    assert!(ScalarEvaluator::from_bytes(family, &corrupt).is_err());
                }
                let other = if family == ScalarIntegral::A0 {
                    ScalarIntegral::B0
                } else {
                    ScalarIntegral::A0
                };
                assert!(ScalarEvaluator::from_bytes(other, &bytes).is_err());
                assert!(
                    restored
                        .evaluate(&[], &mut [Complex::new(0., 0.); 3])
                        .is_err()
                );
                assert!(restored.evaluate_batch(&[], &mut [], usize::MAX).is_err());
                restored.evaluate_batch(&[], &mut [], 0).unwrap();
                let rows: Vec<_> = fixtures.iter().filter(|row| row.family == family).collect();
                let inputs: Vec<_> = rows
                    .iter()
                    .flat_map(|row| row.args.iter().copied())
                    .collect();
                let mut outputs = vec![Complex::new(0., 0.); rows.len() * 3];
                for batch in [1, 3, 4, 5, 31, 256] {
                    outputs.fill(Complex::new(f64::NAN, f64::NAN));
                    for (input, output) in inputs
                        .chunks(family.arity() * batch)
                        .zip(outputs.chunks_mut(3 * batch))
                    {
                        restored
                            .evaluate_batch(input, output, input.len() / family.arity())
                            .unwrap();
                    }
                    for (row, output) in rows.iter().zip(outputs.chunks_exact(3)) {
                        failures.extend(
                            row.failures(output)
                                .into_iter()
                                .map(|failure| format!("batch {batch}: {failure}")),
                        );
                    }
                }
                eprintln!(
                    "{}: {} scalar fixtures, six batch layouts, {} serialized bytes",
                    family.name(),
                    rows.len(),
                    bytes.len()
                );
            }
            assert!(
                failures.is_empty(),
                "{} coefficient comparisons failed:\n{}",
                failures.len(),
                failures.join("\n")
            );
        })
        .unwrap()
        .join()
        .unwrap();
}
