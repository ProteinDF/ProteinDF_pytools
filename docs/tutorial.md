# Tutorial

This tutorial guides you through your first ProteinDF calculation and post-processing with `ProteinDF_pytools`.

## Introduction

In this tutorial, you will:

1. Run a single-point SCF calculation of a water molecule ($H_2O$) with ProteinDF.
2. Archive the calculation results into an HDF5 file (`pdfresults.h5`).
3. Generate plots of total energy, orbital energy levels, and atomic charges.
4. Read calculation results from Python using `PdfParam_H5`.

### Prerequisites

- ProteinDF is installed on your system (providing `PDF.x` and `pdf-setup.sh`).
- `ProteinDF_pytools` installed in your Python environment (see [Installation](installation.md)).

## Setup

### Set environment variables

ProteinDF requires the `PDF_HOME` environment variable pointing to the ProteinDF installation directory. Add the binary directory to your `PATH`:

```bash
export PDF_HOME=/path/to/ProteinDF
export PATH=$PDF_HOME/bin:$PATH
```

### Prepare the working directory

ProteinDF uses a scratch directory named `fl_Work` to store intermediate matrices and vectors during execution. Run `pdf-setup.sh` to create this directory in your working directory:

```bash
pdf-setup.sh
```

### Install pytools

Ensure that `ProteinDF_pytools` is installed in your Python environment according to [Installation](installation.md). After installation, the `pdf-*` command-line tools are available.

## First calculation (H2O)

### Prepare the input file

Create an input file named `fl_Userinput` in your working directory:

```text
>>>>MAIN
	step-control	= [create integral guess scf force]

>>>>INTEGRAL
	cut-value	= 1.0e-10

>>>>GUESS

>>>>SCF
	charge-extrapolate-number	= 0
	scf-start-guess	= core
	orbital_independence_threshold  = 0.00
	orbital_independence_threshold/canonical = 0.00
	orbital_independence_threshold/lowdin = 0.00
	orbital-overlap-correspondence	= off
	orbital-overlap-correspondence-first	= off
	max-iteration	= 300
	method	= rks
	method/rks/electron-number	= 10
	method/rks/occlevel	= [ 1 - 5 ]
	convergence/type	= density
	convergence/threshold	= 1e-3
	convergence/threshold-energy	= 1e-4
	scf-acceleration	= damping
	scf-acceleration/damping/damping-factor	= 0.67
	scf-acceleration/damping/damping-type = density_matrix
	xc-potential	= SVWN
	J_engine = CD
	K_engine = CD
	XC_engine = grid
	CDAM_tau = 1.0E-10
	CD_epsilon = 1.0E-6

>>>>MOLECULE
	geometry/cartesian/unit	= angstrom
	geometry/cartesian/input	= {
		O         0.00            0.00            0.00
		H         1.00            0.00            0.00
		H         0.00            1.00            0.00
	}end

	basis-set/orbital	= {
		H = "O-DZVP2.H"
		O = "O-DZVP2.O"
	}end

	basis-set/density-auxiliary	= {
		H = "A-DZVP2.H"
		O = "A-DZVP2.O"
	}end

	basis-set/exchange-auxiliary	= {
		H = "A-DZVP2.H"
		O = "A-DZVP2.O"
	}end
```

### Convergence criteria

This tutorial explicitly specifies ProteinDF's default convergence criteria:

- `convergence/threshold = 1e-3`
- `convergence/threshold-energy = 1e-4`

These values match the defaults defined in ProteinDF (`data/pdfkwd.yaml`). Using the default criteria allows the calculation for a small test molecule to reach convergence reliably in a few tens of iterations and seconds while demonstrating the standard workflow. (Omitting these two lines in `fl_Userinput` also uses the same default values.)

### Run ProteinDF

Execute `PDF.x` and redirect output to a log file:

```bash
PDF.x > pdf.out 2>&1
```

Check the end of `pdf.out` to confirm convergence:

```text
  21 th TE = -75.8920376510484544
 *** Convergence conditions are satisfied: (1st) ***
  22 th TE = -75.8919893344220782
 *** Convergence conditions are satisfied: (2nd) ***
*** SCF is well converged ***
```

The SCF converged in 22 iterations with a final total energy of `-75.89198933` a.u.

## Archive the results

ProteinDF writes intermediate binary files into `fl_Work`. Use `pdf-archive-h5.py` to collect these files into a single HDF5 archive (`pdfresults.h5`).

Pass the scratch directory with `-w fl_Work`:

```bash
pdf-archive-h5.py -w fl_Work -l full --pop -o pdfresults.h5
```

### Archiving levels

Use `-l <level>` to select which data to store. The default is `standard`. Each higher level includes all data from lower levels:

| Level | Included data |
|---|---|
| `basic` | Molecule, basis sets, control parameters, total energies across all iterations, orbital occupations, overlap matrix S, core Hamiltonian h, dummy atom term h2, converged P, C, and eigenvalues |
| `standard` (default) | `basic` plus input density matrix P_in of iteration 1, energy gradient |
| `full` | `standard` plus P, C, Fock matrix F, and eigenvalues for all iterations |
| `debug` | `full` plus Mulliken charges for all iterations |

- Mulliken charges for the converged solution are saved in `basic`, `standard`, and `full` only when `--pop` is specified. Under `debug`, Mulliken charges for all iterations are saved whether `--pop` is given or not.
- To plot eigenvalue history across iterations, archive with `-l full` or `-l debug` so that eigenvalues of all iterations are stored.

## Plot the results

Use the plotting scripts to visualize the archived results.

### Total energy history

Plot the total energy across SCF iterations:

```bash
pdf-plot-energy.py pdfresults.h5 -o energy.png
```

You can also pass `--yaml` to print the energy values in YAML format.

- *How to read the plot*: The horizontal axis represents the SCF iteration count, and the vertical axis represents the total energy in atomic units (hartree). The energy decreases rapidly in early iterations and flattens out as the SCF reaches convergence.

### Orbital energy level history

Plot how molecular orbital eigenvalues change across iterations:

```bash
pdf-plot-elevel-history.py pdfresults.h5
```

This command outputs `eigvals_hist.rks.png` by default. Note that this requires archiving with `-l full` or `-l debug`.

- *How to read the plot*: The horizontal axis shows iterations, and the vertical axis shows orbital energies in eV. Energy levels are plotted as horizontal lines at each iteration, with the HOMO level highlighted in red and other levels shown in black. You can see orbital levels stabilize as the self-consistent field converges.

### Final iteration energy levels

Generate an inspection report and an energy level diagram for the converged iteration:

```bash
pdf-report-h5.py pdfresults.h5 report
```

This generates `report/eigvals_last.rks.png` and prints a summary of calculation metadata:

```text
=========================
ProteinDF report:
=========================

======================== ================
item                     number
======================== ================
the number of AOs                      24
the number of MOs                      24
the number of iterations               22
Total energy                 -75.89198933
======================== ================
```

You can also plot energy levels directly with `pdf-plot-elevel.py`:

```bash
pdf-plot-elevel.py -d pdfresults.h5 -o elevel.png
```

- *How to read the plot*: Energy levels are drawn vertically. Occupied orbitals are drawn at the bottom and virtual orbitals at the top, showing the HOMO-LUMO gap clearly.

### Mulliken atomic charges

Plot Mulliken atomic charges (when archived with `--pop`):

```bash
pdf-plot-pop.py --h5 pdfresults.h5 -o pop.png
```

- *How to read the plot*: The horizontal axis displays atom indices (1, 2, 3...), and the vertical axis shows the change in Mulliken charge from the initial iteration to the final iteration. Each bar represents the charge difference for that atom.

## Read the results from Python

You can read data directly from `pdfresults.h5` in Python using `PdfParam_H5`.

```python
import proteindf_tools as pdf

# Open the archived HDF5 results
param = pdf.PdfParam_H5("pdfresults.h5")

# Read total energy at the final SCF iteration
last_iter = param.iterations
energy = param.get_total_energy(last_iter)
print(f"Total energy: {energy:.8f} a.u.")

# Read energy gradient matrix (-dE/dR)
# The sign corresponds to force (-dE/dR) in atomic units
grad = param.get_gradient()
print("Gradient matrix (-dE/dR):")
print(grad)
```

Example output:

```text
Total energy: -75.89198933 a.u.
Gradient matrix (-dE/dR):
              1 th       2 th       3 th
  ---- ---------------------------------
     1   -0.035529  -0.035529  -0.000000
     2    0.013237   0.022292   0.000000
     3    0.022292   0.013237   0.000000
```

`get_total_energy(iteration)` returns the total energy in atomic units (hartree) for the specified iteration. `get_gradient()` returns the nuclear gradients as a matrix, where the sign corresponds to $-dE/dR$ (force).

## Troubleshooting

### SCF does not converge

- Increase `max-iteration` in `fl_Userinput`.
- Check if looser convergence thresholds (`convergence/threshold`, `convergence/threshold-energy`) help for initial tests.
- Adjust the SCF damping factor (`scf-acceleration/damping/damping-factor = 0.67`).

### Missing fl_Work directory

- ProteinDF fails if `fl_Work` does not exist before execution. Always run `pdf-setup.sh` before running `PDF.x`.
- When running `pdf-archive-h5.py`, ensure `-w fl_Work` points to the correct scratch directory.

### Cannot open old SQLite .db files

- Older versions of pytools archived results into SQLite databases (`pdfresults.db`).
- Recent tools (`pdf-*-h5.py`, `PdfParam_H5`) use the HDF5 format (`pdfresults.h5`).
- If you have an existing calculation with its `fl_Work` directory intact, re-run `pdf-archive-h5.py` to create a new `pdfresults.h5`.
