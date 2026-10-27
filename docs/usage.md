# Command-line tools

`ProteinDF_pytools` installs a set of `pdf-*` command-line scripts on top
of the `proteindf_tools` library. Each one is a thin wrapper that reads
ProteinDF's binary parameter files (`pdfparam.mpac` / HDF5) and other
result files, and produces a report, plot, or converted file. Run any of
them with `-h` / `--help` for the full list of options.

## Command dispatcher

| Script | Description |
| --- | --- |
| `pdf` (`pdfcmd.py`) | dispatch a `pdf-*` subcommand by name, e.g. `pdf report-h5` runs `pdf-report-h5.py`, `pdf archive-h5` runs `pdf-archive-h5.py` |

## Archiving, reporting & testing

| Script | Description |
| --- | --- |
| `pdf-archive-h5.py` | archive ProteinDF results (HDF5) |
| `pdf-report-h5.py` | make a ProteinDF report (HDF5) |
| `pdf-test.py` | compare ProteinDF results |
| `pdf-test-h5.py` | compare ProteinDF results (HDF5) |
| `pdf-test-eri.py` | test ERI (electron repulsion integral) values |

## Inspecting results

| Script | Description |
| --- | --- |
| `pdf-env.py` | report the ProteinDF runtime environment |
| `pdf-atom-index.py` | output the atom index |
| `pdf-info-geom.py` | output the molecular geometry in XYZ format |
| `pdf-info-iteration.py` | output per-iteration geometry in XYZ format |
| `pdf-info-orb.py` | output orbital information |
| `pdf-info-xyz.py` | output the geometry as an XYZ file |
| `pdf-mat-info.py` | show matrix metadata (size, type, ...) |
| `pdf-mat-show.py` | show matrix contents |
| `pdf-estimate-ao.py` | estimate the AO cell size |
| `pdf-show-profile.py` | [expert] show a profile from a `cProfile` stats file |

## Plotting

| Script | Description |
| --- | --- |
| `pdf-plot-basisset.py` | plot a basis set |
| `pdf-plot-elevel.py` | plot energy levels (single, vertical) |
| `pdf-plot-mo-tracer.py` | plot energy levels across iterations |
| `pdf-plot-mat.py` | plot a matrix |
| `pdf-plot-decaymat.py` | plot a matrix (decay/thinning view) |
| `pdf-plot-dos.py` | plot the density of states |
| `pdf-plot-gvalues.py` | plot G-values |
| `pdf-plot-pop.py` | plot population analysis results |
| `pdf-plot-vector.py` | plot a vector |

## Population analysis & charges

| Script | Description |
| --- | --- |
| `pdf-pop-resp.py` | calculate RESP charges |
| `pdf-pop-classo.py` | fit charges via constrained LASSO regression |
| `pdf-pop-ridge.py` | fit charges via ridge regression |
| `pdf-pop-rrms.py` | calculate the relative RMS (rrms) of a charge fit |
| `pdf-set-charges.py` | write fitted charges into a `pdfparam` file |
| `pdf-reg-harris.py` | register results in the Harris DB |

## Format conversion & input generation

| Script | Description |
| --- | --- |
| `pdf-make-basis2.py` | build a `basis2` basis-set database file |
| `xyz2pdf.py` | make a ProteinDF input from an XYZ file |
| `xyz2gau.py` | make a Gaussian input from an XYZ file |
| `g0xeri2mpac.py` | convert Gaussian output to MessagePack |
