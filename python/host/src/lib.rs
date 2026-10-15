//! Minimal development host: one extension and one Symbolica state.
use pyo3::{prelude::*, types::PyModule};
use symbolica::api::python::{Citation, SymbolicaCommunityModule, create_symbolica_module};

#[pymodule]
fn core(module: &Bound<'_, PyModule>) -> PyResult<()> {
    create_symbolica_module(module)?;
    module.add_function(pyo3::wrap_pyfunction!(get_citations, module)?)?;
    let hep = PyModule::new(module.py(), "symbolica.community.hepkit_native")?;
    oneloop_native::register_hep_module(&hep)?;
    module.add("_hepkit_native", &hep)?;
    let modules = module.py().import("sys")?.getattr("modules")?;
    modules.set_item("symbolica.community.hepkit_native", &hep)?;
    // The legacy import shares the very same classes and Symbolica state.
    modules.set_item(
        "symbolica.community.oneloop_native",
        hep.getattr("_oneloop_native")?,
    )?;
    // Run the community hook; OneLOop numerical backends remain lazy.
    oneloop_native::CommunityModule::initialize(module.py())?;
    Ok(())
}

#[pyo3::pyfunction]
fn get_citations() -> Vec<Citation> {
    let mut citations = vec![Citation {
        id: "doi:10.5281/zenodo.17054381".into(),
        url: "https://doi.org/10.5281/zenodo.17054381".into(),
        reference: "Ben Ruijl. Symbolica (2025). doi:10.5281/zenodo.17054381.".into(),
        bibtex: r#"@software{ruijl_symbolica_2025,
  author = {Ruijl, Ben},
  title = {Symbolica},
  year = {2025},
  version = {0.18.0},
  doi = {10.5281/zenodo.17054381},
  url = {https://zenodo.org/records/17054381}
}"#
        .into(),
        reasons: vec!["Symbolic and numerical computation with Symbolica.".into()],
        description: String::new(),
        relevance: None,
    }];
    citations.extend(oneloop_native::CommunityModule::get_citations());
    citations
}
