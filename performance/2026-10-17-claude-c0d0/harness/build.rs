use std::process::Command;
fn main() {
    let out = std::env::var("OUT_DIR").unwrap();
    let olo = "/shared/external/ninja/references/claude_reference/target/benchmark/oneloop-3.6";
    let st = Command::new("gfortran")
        .args(["-O3", "-fPIC", "-c", "olowrap.f90", "-I", olo, "-J", &out, "-o", &format!("{out}/olowrap.o")])
        .status().unwrap();
    assert!(st.success());
    let st = Command::new("ar").args(["rcs", &format!("{out}/libolowrap.a"), &format!("{out}/olowrap.o")]).status().unwrap();
    assert!(st.success());
    println!("cargo:rustc-link-search=native={out}");
    println!("cargo:rustc-link-lib=static=olowrap");
    println!("cargo:rustc-link-search=native={olo}");
    println!("cargo:rustc-link-lib=static=avh_olo");
    let gf = String::from_utf8(Command::new("gfortran").arg("-print-file-name=libgfortran.so").output().unwrap().stdout).unwrap();
    let dir = std::path::Path::new(gf.trim()).parent().unwrap().to_owned();
    println!("cargo:rustc-link-search=native={}", dir.display());
    println!("cargo:rustc-link-arg=-Wl,-rpath,{}", dir.display());
    println!("cargo:rustc-link-lib=dylib=gfortran");
    println!("cargo:rustc-link-lib=dylib=quadmath");
    println!("cargo:rerun-if-changed=olowrap.f90");
}
