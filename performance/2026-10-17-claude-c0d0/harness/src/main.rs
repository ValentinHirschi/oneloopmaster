//! In-process, interleaved timing of oneloopmaster's native masters and OneLOop 3.6
//! (libavh_olo.a of the o3-fortranrules Ninja build, called through bind(c) wrappers
//! identical to Ninja's avholo_wrapper.F90) on the recorded Ninja master-call streams.
//!
//! usage:
//!   olbench check CASE.masters.json...
//!   olbench time OUT.json CASE.masters.json... [--rounds 9] [--block-ms 25]
//!   olbench percall CASE.masters.json [--reps 200]
use oneloop::{NativeEvaluator, NativeTensorEvaluator, ScalarIntegral, TensorIntegral};
use serde_json::{Value, json};
use std::time::Instant;
use symbolica::domains::float::Complex;

type C = Complex<f64>;

#[repr(C)]
#[derive(Clone, Copy, Default, Debug)]
struct Cd {
    re: f64,
    im: f64,
}

unsafe extern "C" {
    fn olw_d0_cm(r: *mut Cd, p1: *const Cd, p2: *const Cd, p3: *const Cd, p4: *const Cd, p12: *const Cd, p23: *const Cd,
                 m1: *const Cd, m2: *const Cd, m3: *const Cd, m4: *const Cd, rmu: *const f64);
    fn olw_c0_cm(r: *mut Cd, p1: *const Cd, p2: *const Cd, p3: *const Cd, m1: *const Cd, m2: *const Cd, m3: *const Cd, rmu: *const f64);
    fn olw_b0_cm(r: *mut Cd, pp: *const Cd, m1: *const Cd, m2: *const Cd, rmu: *const f64);
    fn olw_b11_cm(b11: *mut Cd, b00: *mut Cd, b1: *mut Cd, b0: *mut Cd, pp: *const Cd, m1: *const Cd, m2: *const Cd, rmu: *const f64);
    fn olw_a0_cm(r: *mut Cd, mm: *const Cd, rmu: *const f64);
    fn olw_setup(thrs: *const f64);
}

#[derive(Clone)]
struct Call {
    // 0 D0, 1 C0, 2 B0, 3 rank-2 bubble (B11 tensor), 4 A0, 5 rank-3 bubble (B111 tensor)
    family: usize,
    input: Vec<C>,
    olo_in: Vec<Cd>,
    rmu: f64,
    reference: Vec<C>,
}

const FAMILY: [&str; 6] = ["D0", "C0", "B0", "B11", "A0", "B111"];
const MAX_OUT: usize = 12;

fn c(v: &Value) -> C {
    C::new(v[0].as_f64().unwrap(), v[1].as_f64().unwrap())
}
fn real(x: f64) -> C {
    C::new(x, 0.0)
}
fn cd(z: C) -> Cd {
    Cd { re: z.re, im: z.im }
}

fn parse_call(v: &Value) -> Option<Call> {
    let kind = v["kind"].as_str().unwrap();
    let s: Vec<f64> = v["s"].as_array().unwrap().iter().map(|x| x.as_f64().unwrap()).collect();
    let m: Vec<C> = v["m2"].as_array().unwrap().iter().map(c).collect();
    let mu2f = v["mu_squared"].as_f64().unwrap();
    let mu2 = real(mu2f);
    let result: Vec<C> = v["result"].as_array().unwrap().iter().map(c).collect();
    let r3 = |o: usize| vec![result[o], result[o + 1], result[o + 2]];
    let rev = |blocks: &[usize]| blocks.iter().flat_map(|&o| r3(o)).collect::<Vec<C>>();
    let (family, input, reference) = match kind {
        "box_cm" => (0, vec![real(s[0]), real(s[1]), real(s[2]), real(s[3]), real(s[4]), real(s[5]), m[0], m[1], m[2], m[3], mu2], r3(0)),
        "triangle_cm" => (1, vec![real(s[0]), real(s[1]), real(s[2]), m[0], m[1], m[2], mu2], r3(0)),
        "bubble_cm" => (2, vec![real(s[0]), m[0], m[1], mu2], r3(0)),
        "rank2_bubble_cm" => (3, vec![real(s[0]), m[0], m[1], mu2], rev(&[6, 3, 0])),
        "rank3_bubble_cm" => (5, vec![real(s[0]), m[0], m[1], mu2], rev(&[9, 6, 3, 0])),
        "tadpole_cm" => (4, vec![m[0], mu2], r3(0)),
        _ => return None,
    };
    let olo_in = input[..input.len() - 1].iter().map(|&z| cd(z)).collect();
    Some(Call { family, input, olo_in, rmu: mu2f.sqrt(), reference })
}

struct Native {
    ev: Vec<NativeEvaluator<f64>>,
    b11: NativeTensorEvaluator<f64>,
    b111: NativeTensorEvaluator<f64>,
}

impl Native {
    fn new() -> Self {
        let fam = [ScalarIntegral::D0, ScalarIntegral::C0, ScalarIntegral::B0, ScalarIntegral::B0, ScalarIntegral::A0];
        Self {
            ev: fam.iter().map(|&f| NativeEvaluator::<f64>::new(f).unwrap()).collect(),
            b11: NativeTensorEvaluator::<f64>::new(TensorIntegral::B11).unwrap(),
            b111: NativeTensorEvaluator::<f64>::new(TensorIntegral::B111).unwrap(),
        }
    }
    #[inline]
    fn eval(&mut self, call: &Call, out: &mut [C; MAX_OUT]) {
        let k = call.reference.len();
        match call.family {
            3 => self.b11.evaluate(&call.input, &mut out[..k]).unwrap(),
            5 => self.b111.evaluate(&call.input, &mut out[..k]).unwrap(),
            f => self.ev[f].evaluate(&call.input, &mut out[..k]).unwrap(),
        }
    }
}

/// OneLOop through the same entry points as Ninja's AvHOneLoop (rank-3 bubbles as
/// IntegralLibrary::getRank3BubbleIntegralCM). Output in the native order.
#[inline]
fn olo_eval(call: &Call, out: &mut [C; MAX_OUT]) {
    let x = &call.olo_in;
    let mut r = [Cd::default(); 12];
    unsafe {
        match call.family {
            0 => olw_d0_cm(r.as_mut_ptr(), &x[0], &x[1], &x[2], &x[3], &x[4], &x[5], &x[6], &x[7], &x[8], &x[9], &call.rmu),
            1 => olw_c0_cm(r.as_mut_ptr(), &x[0], &x[1], &x[2], &x[3], &x[4], &x[5], &call.rmu),
            2 => olw_b0_cm(r.as_mut_ptr(), &x[0], &x[1], &x[2], &call.rmu),
            4 => olw_a0_cm(r.as_mut_ptr(), &x[0], &call.rmu),
            3 | 5 => {
                let (mut b11, mut b00, mut b1, mut b0) = ([Cd::default(); 3], [Cd::default(); 3], [Cd::default(); 3], [Cd::default(); 3]);
                olw_b11_cm(b11.as_mut_ptr(), b00.as_mut_ptr(), b1.as_mut_ptr(), b0.as_mut_ptr(), &x[0], &x[1], &x[2], &call.rmu);
                let z = |a: Cd| C::new(a.re, a.im);
                for k in 0..3 {
                    out[k] = z(b0[k]);
                    out[3 + k] = z(b1[k]);
                    out[6 + k] = z(b11[k]);
                }
                if call.family == 5 {
                    rank3(call, &b11, &b1, &b0, out);
                }
                return;
            }
            _ => unreachable!(),
        }
    }
    for k in 0..3 {
        out[k] = C::new(r[k].re, r[k].im);
    }
}

#[inline(never)]
fn rank3(call: &Call, b11: &[Cd; 3], b1: &[Cd; 3], b0: &[Cd; 3], out: &mut [C; MAX_OUT]) {
    let z = |a: Cd| C::new(a.re, a.im);
    let x = &call.olo_in;
    let p2 = x[0].re;
    let (m1q, m2q) = (z(x[1]), z(x[2]));
    let rl = |v: f64| C::new(v, 0.0);
    let mut b111 = [C::new(0.0, 0.0); 3];
    if p2.abs() > 1e-8 {
        let mut a1 = [Cd::default(); 3];
        let mut a2 = [Cd::default(); 3];
        unsafe {
            olw_a0_cm(a1.as_mut_ptr(), &x[1], &call.rmu);
            olw_a0_cm(a2.as_mut_ptr(), &x[2], &call.rmu);
        }
        let mut i2 = [C::new(0.0, 0.0); 3];
        i2[0] = z(a2[0]) + m1q * z(b0[0]) - (rl(p2) - rl(3.0) * (m1q + m2q)) / rl(6.0);
        i2[1] = z(a2[1]) + m1q * z(b0[1]);
        i2[2] = z(a2[2]) + m1q * z(b0[2]);
        let b00: Vec<C> = (0..3).map(|k| (i2[k] - rl(p2) * z(b11[k])) * rl(0.25)).collect();
        let rat = -rl(p2 * p2) * (rl(p2) - rl(2.0) * m1q - rl(4.0) * m2q) / rl(3.0) + rl(0.5 * p2) * (m1q * m1q - m2q * m2q);
        let c1 = -rl(4.0 * p2 * p2) * m1q;
        let c2 = rl(4.0) * (m2q - m1q - rl(p2));
        for k in 0..3 {
            b111[k] = c1 * z(b1[k]) + c2 * rl(p2) * (b00[k] + rl(p2) * z(b11[k])) + rl(p2) * (m1q * z(a1[k]) - m2q * z(a2[k]));
            if k == 0 {
                b111[k] += rat;
            }
            b111[k] = b111[k] / rl(4.0 * p2 * p2 * p2);
        }
    } else {
        let mut z11 = [Cd::default(); 3];
        let zero = Cd::default();
        unsafe { olw_b0_cm(z11.as_mut_ptr(), &zero, &x[1], &x[1], &call.rmu) };
        if m1q != m2q {
            let mut z22 = [Cd::default(); 3];
            unsafe { olw_b0_cm(z22.as_mut_ptr(), &zero, &x[2], &x[2], &call.rmu) };
            let a3 = m1q - m2q;
            let a3q = a3 * a3;
            let a3c = a3q * a3;
            let a4 = rl(1.0) / (a3c * a3);
            let m1qq = m1q * m1q;
            let m1qc = m1qq * m1q;
            let m1qf = m1qq * m1qq;
            let m2qq = m2q * m2q;
            let m2qc = m2qq * m2q;
            let bb1 = -(m1qf * a4) / rl(4.0);
            let bb2 = (rl(-1.0) + m1qf * a4) / rl(4.0);
            let bb3 = (rl(-25.0) * m1qc + rl(23.0) * m1qq * m2q - rl(13.0) * m1q * m2qq + rl(3.0) * m2qc) / (rl(48.0) * a3c);
            for k in 0..3 {
                b111[k] = bb1 * z(z11[k]) + bb2 * z(z22[k]);
            }
            b111[0] += bb3;
        } else {
            for k in 0..3 {
                b111[k] = rl(-0.25) * z(z11[k]);
            }
        }
    }
    for k in 0..3 {
        out[9 + k] = b111[k];
    }
}

fn median(mut v: Vec<f64>) -> f64 {
    v.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let n = v.len();
    if n % 2 == 1 { v[n / 2] } else { 0.5 * (v[n / 2 - 1] + v[n / 2]) }
}

fn load(file: &str) -> (String, Vec<Vec<Call>>) {
    let data: Value = serde_json::from_str(&std::fs::read_to_string(file).unwrap()).unwrap();
    let lanes = data["lanes"].as_array().unwrap().iter()
        .map(|lane| lane["calls"].as_array().unwrap().iter().filter_map(parse_call).collect())
        .collect();
    (data["case"].as_str().unwrap().to_string(), lanes)
}

fn sdiff(a: &[C], b: &[C]) -> f64 {
    let scale = b.iter().map(|z| z.re.hypot(z.im)).fold(0.0, f64::max).max(1e-300);
    a.iter().zip(b).map(|(x, y)| (x.re - y.re).hypot(x.im - y.im) / scale).fold(0.0, f64::max)
}

/// The 8 dihedral relabellings of a box / 6 of a triangle (input order as NativeEvaluator).
fn relabellings(call: &Call) -> Vec<Vec<C>> {
    let x = &call.input;
    match call.family {
        1 => {
            // p1=(12), p2=(23), p3=(31); masses m1,m2,m3. permutations of vertices.
            let perms = [[0, 1, 2], [1, 2, 0], [2, 0, 1], [0, 2, 1], [2, 1, 0], [1, 0, 2]];
            perms.iter().map(|q| {
                let m = [x[3], x[4], x[5]];
                // p between internal lines i and j: p(0,1)=x0, p(1,2)=x1, p(2,0)=x2
                let p = |i: usize, j: usize| -> C {
                    let (a, b) = (i.min(j), i.max(j));
                    match (a, b) { (0, 1) => x[0], (1, 2) => x[1], (0, 2) => x[2], _ => unreachable!() }
                };
                vec![p(q[0], q[1]), p(q[1], q[2]), p(q[2], q[0]), m[q[0]], m[q[1]], m[q[2]], x[6]]
            }).collect()
        }
        0 => {
            let perms = [[0, 1, 2, 3], [1, 2, 3, 0], [2, 3, 0, 1], [3, 0, 1, 2], [0, 3, 2, 1], [3, 2, 1, 0], [2, 1, 0, 3], [1, 0, 3, 2]];
            perms.iter().map(|q| {
                let m = [x[6], x[7], x[8], x[9]];
                let p = |i: usize, j: usize| -> C {
                    let (a, b) = (i.min(j), i.max(j));
                    match (a, b) { (0, 1) => x[0], (1, 2) => x[1], (2, 3) => x[2], (0, 3) => x[3], (0, 2) => x[4], (1, 3) => x[5], _ => unreachable!() }
                };
                vec![p(q[0], q[1]), p(q[1], q[2]), p(q[2], q[3]), p(q[3], q[0]), p(q[0], q[2]), p(q[1], q[3]),
                     m[q[0]], m[q[1]], m[q[2]], m[q[3]], x[10]]
            }).collect()
        }
        _ => vec![],
    }
}

fn check(files: &[String]) {
    let mut nat = Native::new();
    let mut out = [C::new(0.0, 0.0); MAX_OUT];
    let mut alt = [C::new(0.0, 0.0); MAX_OUT];
    let mut olo = [C::new(0.0, 0.0); MAX_OUT];
    let mut worst_all = [0.0f64; 6];
    for file in files {
        let (name, lanes) = load(file);
        let mut md = [0.0f64; 6];
        let mut mo = [0.0f64; 6];
        let mut ms = [0.0f64; 6];
        let mut worst: Vec<Option<String>> = vec![None; 6];
        for lane in &lanes {
            for call in lane {
                let k = call.reference.len();
                nat.eval(call, &mut out);
                olo_eval(call, &mut olo);
                let d = sdiff(&out[..k], &call.reference);
                let dol = sdiff(&olo[..k], &call.reference);
                if dol > mo[call.family] { mo[call.family] = dol; }
                if d > md[call.family] {
                    md[call.family] = d;
                    worst[call.family] = Some(format!("in={:?}\n      nat={:?}\n      olo={:?}",
                        call.input.iter().map(|z| (z.re, z.im)).collect::<Vec<_>>(),
                        out[..k].iter().map(|z| (z.re, z.im)).collect::<Vec<_>>(),
                        call.reference.iter().map(|z| (z.re, z.im)).collect::<Vec<_>>()));
                }
                for perm in relabellings(call) {
                    let c2 = Call { input: perm, ..call.clone() };
                    nat.eval(&c2, &mut alt);
                    let d = sdiff(&alt[..3], &out[..3]);
                    if d > ms[call.family] { ms[call.family] = d; }
                }
            }
        }
        print!("{name:<28}");
        for f in 0..6 {
            if md[f] > 0.0 || ms[f] > 0.0 {
                print!(" {}:{:.1e}/rl{:.0e}", FAMILY[f], md[f], ms[f]);
            }
            worst_all[f] = worst_all[f].max(md[f]);
        }
        println!("   (live OneLOop vs recorded: {:.1e})", mo.iter().cloned().fold(0.0, f64::max));
        if std::env::var("WORST").is_ok() {
            for f in 0..2 {
                if let Some(w) = &worst[f] { println!("   worst {}: {}", FAMILY[f], w); }
            }
        }
    }
    println!("max per family: {:?}", FAMILY.iter().zip(worst_all.iter()).map(|(a, b)| format!("{a}:{b:.1e}")).collect::<Vec<_>>());
}

fn time_group(sel: &[Vec<Call>], nat: &mut Native, block_ms: f64, rounds: usize) -> (f64, f64, Vec<f64>, Vec<f64>) {
    let mut out = [C::new(0.0, 0.0); MAX_OUT];
    let mut sink = 0.0f64;
    let mut run_nat = |iters: u64, nat: &mut Native| -> f64 {
        let t0 = Instant::now();
        for _ in 0..iters {
            for lane in sel {
                for call in lane {
                    nat.eval(call, &mut out);
                    sink += out[0].re;
                }
            }
        }
        t0.elapsed().as_nanos() as f64
    };
    let mut out2 = [C::new(0.0, 0.0); MAX_OUT];
    let mut sink2 = 0.0f64;
    let mut run_olo = |iters: u64| -> f64 {
        let t0 = Instant::now();
        for _ in 0..iters {
            for lane in sel {
                for call in lane {
                    olo_eval(call, &mut out2);
                    sink2 += out2[0].re;
                }
            }
        }
        t0.elapsed().as_nanos() as f64
    };
    let calib = |f: &mut dyn FnMut(u64) -> f64| -> u64 {
        let mut iters = 1u64;
        loop {
            let t = f(iters);
            if t >= 0.2 * block_ms * 1e6 {
                return ((block_ms * 1e6) / (t / iters as f64)).ceil().max(1.0) as u64;
            }
            iters *= 4;
        }
    };
    let it_n = calib(&mut |i| run_nat(i, nat));
    let it_o = calib(&mut |i| run_olo(i));
    run_nat(it_n, nat);
    run_olo(it_o);
    let pn = it_n as f64 * sel.len() as f64;
    let po = it_o as f64 * sel.len() as f64;
    let mut bn = vec![];
    let mut bo = vec![];
    for _ in 0..rounds {
        bn.push(run_nat(it_n, nat) / pn);
        bo.push(run_olo(it_o) / po);
    }
    std::hint::black_box((sink, sink2));
    (median(bn.clone()), median(bo.clone()), bn, bo)
}

fn time(out_file: &str, files: &[String], block_ms: f64, rounds: usize, only: Option<String>) {
    let mut nat = Native::new();
    let mut reports = vec![];
    let groups: [(&str, &dyn Fn(&Call) -> bool); 5] = [
        ("all", &|_c: &Call| true),
        ("D0", &|c: &Call| c.family == 0),
        ("C0", &|c: &Call| c.family == 1),
        ("B", &|c: &Call| c.family == 2 || c.family == 3 || c.family == 5),
        ("A0", &|c: &Call| c.family == 4),
    ];
    println!("{:<28} {:>9} {:>9} {:>6} | {:>8} {:>8} {:>6} | {:>8} {:>8} {:>6} | {:>7} {:>7} | {:>6} {:>6}",
             "case (us/point)", "OL all", "nat all", "ratio", "OL D0", "nat D0", "ratio", "OL C0", "nat C0", "ratio", "OL B", "nat B", "OL A0", "nat A0");
    for file in files {
        let (name, lanes) = load(file);
        let mut timing = serde_json::Map::new();
        let mut row = vec![];
        for (g, f) in groups.iter() {
            if let Some(o) = &only {
                if g != o && *g != "all" { row.push(None); continue; }
            }
            let sel: Vec<Vec<Call>> = lanes.iter().map(|l| l.iter().filter(|c| f(c)).cloned().collect()).collect();
            if sel.iter().all(|l| l.is_empty()) {
                row.push(None);
                continue;
            }
            let (n, o, bn, bo) = time_group(&sel, &mut nat, block_ms, rounds);
            timing.insert(g.to_string(), json!({"native_ns": n, "oneloop_ns": o, "ratio": n / o, "native_blocks": bn, "oneloop_blocks": bo}));
            row.push(Some((n, o)));
        }
        let f = |x: &Option<(f64, f64)>, i: usize| x.map(|v| format!("{:.2}", if i == 0 { v.1 } else { v.0 } / 1e3)).unwrap_or("-".into());
        let r = |x: &Option<(f64, f64)>| x.map(|v| format!("{:.2}x", v.0 / v.1)).unwrap_or("-".into());
        println!("{:<28} {:>9} {:>9} {:>6} | {:>8} {:>8} {:>6} | {:>8} {:>8} {:>6} | {:>7} {:>7} | {:>6} {:>6}",
                 name, f(&row[0], 0), f(&row[0], 1), r(&row[0]), f(&row[1], 0), f(&row[1], 1), r(&row[1]),
                 f(&row[2], 0), f(&row[2], 1), r(&row[2]), f(&row[3], 0), f(&row[3], 1), f(&row[4], 0), f(&row[4], 1));
        reports.push(json!({"case": name, "file": file, "timing_ns_per_point": timing}));
    }
    let out = json!({"format": "olbench-v1", "block_ms": block_ms, "rounds": rounds, "cases": reports});
    std::fs::write(out_file, serde_json::to_string_pretty(&out).unwrap()).unwrap();
}

fn percall(file: &str, reps: usize, family: Option<usize>) {
    let mut nat = Native::new();
    let (_, lanes) = load(file);
    let mut out = [C::new(0.0, 0.0); MAX_OUT];
    for call in &lanes[0] {
        if call.family > 1 { continue; }
        if let Some(f) = family { if call.family != f { continue; } }
        let mut tn = vec![];
        let mut to = vec![];
        for _ in 0..7 {
            let t0 = Instant::now();
            for _ in 0..reps { nat.eval(call, &mut out); std::hint::black_box(&out); }
            tn.push(t0.elapsed().as_nanos() as f64 / reps as f64);
            let t0 = Instant::now();
            for _ in 0..reps { olo_eval(call, &mut out); std::hint::black_box(&out); }
            to.push(t0.elapsed().as_nanos() as f64 / reps as f64);
        }
        let (n, o) = (median(tn), median(to));
        println!("{} nat {:>8.0} ns  olo {:>8.0} ns  {:.2}x  s={:?} m={:?}", FAMILY[call.family], n, o, n / o,
                 call.input[..if call.family == 0 { 6 } else { 3 }].iter().map(|z| z.re).collect::<Vec<_>>(),
                 call.input[if call.family == 0 { 6..10 } else { 3..6 }].iter().map(|z| (z.re, z.im)).collect::<Vec<_>>());
    }
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    oneloop::initialize().unwrap();
    unsafe { olw_setup(&1.0e-9) };
    let mut block_ms = 25.0;
    let mut rounds = 9usize;
    let mut reps = 200usize;
    let mut only = None;
    let mut family = None;
    let mut files = vec![];
    let mut i = 2;
    while i < args.len() {
        match args[i].as_str() {
            "--block-ms" => { block_ms = args[i + 1].parse().unwrap(); i += 1; }
            "--rounds" => { rounds = args[i + 1].parse().unwrap(); i += 1; }
            "--reps" => { reps = args[i + 1].parse().unwrap(); i += 1; }
            "--only" => { only = Some(args[i + 1].clone()); i += 1; }
            "--family" => { family = Some(args[i + 1].parse().unwrap()); i += 1; }
            f => files.push(f.to_string()),
        }
        i += 1;
    }
    match args[1].as_str() {
        "check" => check(&files),
        "time" => time(&files[0], &files[1..], block_ms, rounds, only),
        "percall" => percall(&files[0], reps, family),
        m => panic!("unknown mode {m}"),
    }
}
