# ProteinDF_pytools

`ProteinDF_pytools` (package name `proteindf_tools`) is a Python library
and set of command-line tools for pre- and post-processing around the
[ProteinDF](https://proteindf.github.io/) quantum chemistry engine.

It provides:

- Python classes to read and manipulate ProteinDF's binary parameter files
  (`pdfparam.mpac` / HDF5), matrices, vectors, basis sets, and orbital
  information (`PdfParam`, `Matrix`, `BasisSet`, `OrbInfo`, ...).
- A collection of `pdf-*` command-line tools to inspect, plot, and report
  on ProteinDF calculation results without writing any Python code.

```{toctree}
:maxdepth: 2
:caption: Contents

installation
usage
api/modules
```

## Indices and tables

- {ref}`genindex`
- {ref}`modindex`
- {ref}`search`
