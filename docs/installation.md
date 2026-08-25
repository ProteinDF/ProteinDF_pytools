# Installation

## Requirements

- Python 3.7 or later
- [`ProteinDF_bridge`](https://proteindf.github.io/ProteinDF_bridge/)
- `numpy`, `scipy`, `pandas`, `matplotlib`, `scikit-learn`, `h5py`,
  `BeautifulSoup4`, `rich`

## From source

```bash
git clone https://github.com/ProteinDF/ProteinDF_pytools.git
cd ProteinDF_pytools
python -m venv .venv
source .venv/bin/activate
pip install .
```

This also installs the `pdf-*` command-line tools listed in
[Command-line tools](usage.md).
