# Consolidated main validation

Validation started on 2026-10-16 with Rust 1.98.1, Linux x86_64,
`RUST_MIN_STACK=134217728`, at most eight Cargo jobs, and serial test execution.
The consolidated implementation is based on `1d27666`, including the existing
PR merge ancestry and native batch emitter commit `79e5f8b`.

The repository owner explicitly requested direct updates to oneloopmaster's
`main` and consolidation of its existing pull requests, without opening new
oneloopmaster pull requests. The consolidation retains the heads and merge
ancestry of PRs 3, 4, 5, 6, 7 and 8. This instruction does not change other
upstream repositories' contribution policies. The Symbolica source pin below
publishes an unchanged existing integration; it introduces no new Symbolica
implementation or pull request.

## One shared Symbolica source

All three standalone build roots (`Cargo.toml`, `python/Cargo.toml`, and
`python/host/Cargo.toml`) patch Symbolica, Numerica and Graphica to
`8e23925ed7b5cbb91607c8633dda114cb43be8ce` from
<https://github.com/ValentinHirschi/symbolica>. This integration contains the
recursive shared-function Dual support, vector branch-scoping repair and
polylog-at-zero series needed by the coincident-dilogarithm constructor.
The common package version `3.0.1` does not by itself identify these fixes.

`cargo metadata --format-version 1`, repeated with each Python manifest,
resolved exactly one copy of each of the three packages, all from that exact
git revision. The fetched git checkout and the clean local checkout used for
the initial Rust test build have identical git tree
`421cd3458040550effaa861d8382267242968c2f`. Thus the local path build and the
published git pin use identical source content, though Cargo assigns different
source identities to their artifacts.

Both tracked Python lockfiles were refreshed. The initial changes replace the three
registry package sources with this git pin, and update Graphica from 3.0.0 to
the integration's 3.0.1. No other package versions changed. The following
standalone checks used the actual git patches, without local source overrides:

```sh
cargo check --locked -j 4 --manifest-path python/Cargo.toml
cargo check --locked -j 4 --manifest-path python/host/Cargo.toml
```

Both passed. These are compilation checks; they do not claim to test Python
extension loading or import behavior.

## Rust tests and retained limitations

The initial compile used all features and built all 35 unit/integration test
executables successfully. Their runtime checks are executed serially, with one
process per target, to preserve each test target's initialization semantics.
The corresponding ordinary standalone commands are:

```sh
export RUST_MIN_STACK=134217728
cargo test -j 8 --all-features --lib --tests --no-run
cargo test -j 8 --all-features --lib --tests -- --nocapture --test-threads=1
```

The unit target passed 78 tests, with one pre-existing ignored construction-cost
diagnostic. Focused integration targets also passed: eight function-definition
tests, 11 native-primitive tests, and the compiled native batch emitter probe.
Together with the 12 emitter and three native workspace/constants unit tests,
this reruns all 27 checks from the original batch-emitter validation against
the consolidated source and the newer Symbolica integration. The generated-code
probe retained its instruction-major trace, zero warmed f64 allocation checks,
Float precision preservation and higher/mixed flattened-Dual comparisons.
It completed in 24.03 seconds, including its generated Rust compilation.
The full-expansion target passed all three tests. Formatting and
`git diff --check` passed.

The complete pre-fix run finished all 35 targets: 169 tests passed, six failed,
and six were ignored. All six failures reported the same D0 JIT argument-limit
error described below. The ignored tests include fresh-process helpers run by
their parent tests, as well as existing diagnostic tests; they are not six
additional passing checks. The original batch-emitter probe also passed when
reached by this complete serial run.

The `evaluators`, `initialization`, `master_callbacks_jit`, `native_hook_errors`,
`native_hooks`, and `native_startup` targets reproduced the eager D0
initialization failure:

```text
could not initialize D0 evaluator: Could not JIT-compile sub-evaluator
'oneloop___olo_boxc_box': function call with 70 args, but only 64 are
available. Use `Config.set_num_args`.
```

Before this failure, the same target passed its A0 (4), B0 (13), dB0 (12), and
C0 (125) scalar fixture comparisons and six batch layouts for each family.
The D0 failure was previously reproduced independently on the unmodified
scalar-emitter base `de0114a`, as recorded in
[native batch emitter validation](native-batch-emitter-validation.md).
The pre-fix failures are retained above rather than counted as passing tests.

## JIT argument-storage follow-up

Commit `384b6375fc01fe5de03fd4cd1116f22803e3dffe`, merged as
`1bdcb2893db667744c437e37add4b7212ff8ab60`, derives the required JIT argument
storage from the existing exported instructions, including nested complex
helper bodies. It uses SymJIT's public configuration API and raises the compatible
minimum SymJIT version to 2.27.3. No integral expressions or numerical tolerances
change.

The repair's two new unit tests passed, as did the eager-initialization and
portable-evaluator tests. The latter covers all 293 original scalar fixtures,
including 139 D0 cases, six batch layouts and serialized-cache restoration.
The initialization target passed its parent test and all six fresh-process
entry routes. The bare-master JIT target passed all 293 fixtures in five batch
layouts, six Laurent-tag orders and original/restored JIT instances.

The three remaining affected targets also passed after the merge in the
consolidated worktree: `native_hook_errors` (one test), `native_hooks` (one test),
and `native_startup` (one parent test covering all three fresh-process entry
routes, plus its intentionally ignored child helper). Thus every failure from
the complete pre-fix run has a passing focused rerun after the repair. The
unaffected targets were not redundantly rerun; this combined evidence is not
reported as a new single full-suite execution.

Both Python `cargo check --locked -j 4` commands passed again after the repair,
using the actual git-pinned Symbolica integration. Their metadata again resolves
exactly one Symbolica, Numerica and Graphica, and now one SymJIT 2.27.3. The
follow-up Python lockfile changes update SymJIT from 2.26.0 to 2.27.3 and add its
`corosensei`/`scopeguard` dependencies.

The existing formatting-only wrap in `src/epsilon.rs` changes no expression or
algorithm. It is included so repository-wide formatting checks are reproducible.

## Evidence files

The validation host retains command logs and per-target return codes under
`/tmp/ninjars-oneloop-main-validation/`: `dependency-resolution.json`,
`root-metadata.json`, `python-metadata.json`, `host-metadata.json`,
`python-check.log`, `host-check.log`, `results.json`, and individual test logs.
The focused targets have separate `priority-results.json` and `priority-*.log`
records.
The JIT follow-up uses `jit-followup-results.json`, `jit-followup-*.log`,
`followup-dependency-resolution.json`, and both `*-followup-check.log` files.
The repair's initial focused logs are `/tmp/ninjars-jit-arity-evaluators.log`
and `/tmp/ninjars-jit-arity-additional.log`.
The build log is `/tmp/ninjars-oneloop-main-build.log`. These temporary host
paths are audit evidence, not prerequisites for the commands above.
