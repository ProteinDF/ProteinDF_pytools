# ProteinDF tools -- pre/post scripts for the ProteinDF

## How to install 'ProteinDF tools'

### notice

The package name has been changed from 2017.04. Please be careful.


### install

1. Clone a copy of the main ProteinDF tools git repo by running:

```bash
git clone git://github.com/ProteinDF/ProteinDF_pytools.git
```

1. Enter the source directory and run build script:

```bash
$ cd ProteinDF_pytools
$ pip install --use-feature=in-tree-build .
```


## Documents

Documentation (English / 日本語) is published at
https://proteindf.github.io/ProteinDF_pytools/ via GitHub Pages.

It can also be built locally:

```bash
pip install -e .[docs]
make -C docs html      # English HTML -> docs/_build/html/en
make -C docs html-ja   # Japanese HTML -> docs/_build/html/ja
make -C docs latexpdf  # English PDF (requires a local TeX install)
```


## License

'ProteinDF tools' is licensed under the GNU GPL v3.
The source code can be found on the Github.


## Bugs

If you find any bugs, please let me know.
And if you have a suggestion for improvement, please let me know.
