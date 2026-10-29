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

## HDF5 に保存されるデータ(`pdf-archive-h5.py`)

### 保存レベル

`pdf-archive-h5.py -l <レベル>` で、保存する範囲を選ぶ。既定は `standard`。各レベルは、前のレベルを含む。

| レベル | 含むもの |
|---|---|
| `basic` | 分子・基底・制御、全反復の全エネルギー、占有数、重なり行列 S、コア Hamiltonian h、ダミー原子の項 h2、収束解の P・C・固有値 |
| `standard`(既定) | `basic` に加えて、反復 1 の入力密度 P_in、エネルギー勾配 gradient |
| `full` | `standard` に加えて、全反復の P・C・F・固有値 |
| `debug` | `full` に加えて、全反復の Mulliken 電荷 |

- Mulliken 電荷は、`--pop` を付けたときだけ、収束解 1 本を保存します(`basic`・`standard`・`full`)。`debug` では、`--pop` の有無にかかわらず、全反復を保存します。

### 保存の項目

| 項目 | 内容 | キー | 含まれるレベル |
|---|---|---|---|
| 分子・基底・制御 | 計算の設定 | `molecule`、`basisset`、`control` | 全レベル |
| 全エネルギー | 反復ごとの全エネルギー(1 本に全反復分) | `TEs` | 全レベル |
| 占有数 | 軌道の占有数(run type ごと) | `occ/<run_type>` | 全レベル |
| 重なり行列 | S | `s` | 全レベル |
| コア Hamiltonian | h(Hpq、一電子のコア項) | `h` | 全レベル |
| ダミー原子の項 | h2(Hpq2、ダミー原子の点電荷による一電子項。ダミー原子がない系では 0) | `h2` | 全レベル |
| 出力密度 P | SCF の密度行列 | `P/<run_type>_<N>` | `basic`・`standard` は収束解、`full`・`debug` は全反復 |
| C 行列 | 分子軌道の係数 | `C/<run_type>_<N>` | `basic`・`standard` は収束解、`full`・`debug` は全反復 |
| 固有値 | 軌道のエネルギー準位 | `energy_level/<run_type>_<N>` | `basic`・`standard` は収束解、`full`・`debug` は全反復 |
| 入力密度 P_in | SCF の反復 1 の入力密度 | `P_in/<run_type>_1` | `standard` 以上 |
| エネルギー勾配 | dE/dR(原子数 × 3)。ProteinDF の `force` の符号を反転したもの | `gradient` | `standard` 以上 |
| F 行列 | Fock 行列 | `F/<run_type>_<N>` | `full`・`debug`(全反復) |
| Mulliken 電荷 | 原子ごとの電荷 | `pop/mulliken/<run_type>_<N>` | `--pop` で収束解 1 本。`debug` は全反復 |

- `<N>` は反復番号です(1 から始まります)。
- 収束解は、最後の反復(`iterations`)です。

### run type の対応

| 方式 | 密度(P・P_in・占有数・Mulliken) | 軌道(固有値・C・F) |
|---|---|---|
| RKS | `rks` | `rks` |
| UKS | `uks_alpha`、`uks_beta` | `uks_alpha`、`uks_beta` |
| ROKS | `roks_closed`、`roks_open` | `roks` |
