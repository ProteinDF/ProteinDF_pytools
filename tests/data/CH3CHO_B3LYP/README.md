# CH3CHO / B3LYP テストデータ

アセトアルデヒド(CH3CHO)の B3LYP 1 点計算の入力と結果。

| ファイル | 内容 |
|---|---|
| `fl_Userinput` | ProteinDF の入力(B3LYP、グリッド 7 相当)。 |
| `reference.h5` | 基準の結果(2026-10-05)。PySCF との全エネルギーの差は -2.15e-5 Ha。**基準はこのファイル。** |
| `pdfresults.h5` | ProteinDF の結果(2026-10-17、反復 51)。`reference.h5` との差は 8.6e-6 Ha。 |

## 入れなかったもの

- `fl_Work/`(ProteinDF の中間ファイル、約 600 件): 再計算で作れるため。
- `pdfresults_std.db`(SQLite の古い結果): SQLite 版は廃止済み。
- `fl_Out_Std`、`pdfparam.mpac`、`force.mtx`、`pdf-archive-h5.log`: 再計算や `pdf-archive-h5.py` で作れるため。
- `.envrc`: この環境の PDF_HOME の設定で、個人のパスを含むため。
