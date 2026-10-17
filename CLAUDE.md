# Claude(レビュー担当)向けルール

このリポジトリでは、**Claudeはレビュー担当**であり、実装はagy / GitHub Copilotに任せる。全体の流れは`CONTRIBUTING.md`、実装担当のルールは`AGENTS.md`を参照。

## Claudeがやること

1. **タスク指示書を書く**: ユーザーと相談して`private/tasks/TASK_<name>.md`を`devtool/TASK_TEMPLATE.md`から作る。**タスク指示書・`TODO.md`・`SPEC.md`は非公開**で、ProteinDF_pytools本体とは別の非公開リポジトリで管理する(developのworktreeの`private/`はそこへのシンボリックリンク、`TODO.md`・`SPEC.md`も`private/`へのリンク)。書いたら非公開リポジトリにコミットする(`git -C private commit`)。ProteinDF_pytools本体にはコミットしない。不具合の分析などを含むため、内容をコミットメッセージや公開されるファイルに書き写さない。設計判断が必要な点は書く前にユーザーに確認する。
   - **指示書は短く、タスクは小さくする**: 指示書は約3KB、対象・完了の定義とも各5項目以内を目安にし、背景は5行以内、調査の経緯や数値は`TODO.md`・`SPEC.md`に置く(実装担当は毎回これを読み込むのでコンテキストを使う)。超えるなら調査・修正・基準値の更新などに分けて複数のTASKにする。レビュー結果も要点だけを書く。
2. **実装を依頼する**: ユーザーの指示があれば`devtool/delegate.sh private/tasks/TASK_<name>.md`を実行する(時間がかかるのでバックグラウンドで実行する)。ユーザー自身がagy / Copilotに渡す場合もある。
   - agyが利用上限で止まると`delegate.sh`は終了コード75で終わる。再開の方法(上限のリセットを待って`-c`で続ける、`-W`で自動的に待つ、`-a copilot`に切り替える)はユーザーに確認する。
3. **レビューする**: ブランチのworktree(`devtool/flow.sh path <branch>`)で`/code-review`を使い、あわせて次を行う。
   - `devtool/review.sh <branch>`を**自分で実行する**(`check.sh`を実行し、出力を`.git/review-logs/`に保存する)。実装担当の完了報告の数字・出力を鵜呑みにしない。
   - タスク指示書の「完了の定義」を1項目ずつ照合する。
   - スコープ外の変更、警告の抑制、エラーの握りつぶし、変更していない行の整形がないか確認する。
4. **レビュー結果を書く**: タスク指示書の末尾に`## レビュー結果(N回目、YYYY-MM-DD、要修正|収束)`として追記し、非公開リポジトリにコミットする。要修正なら番号付きの修正依頼と「完了の定義(修正後)」を書く。実装担当は同じブランチに追加コミットする。
5. **マージする**: 収束したら、**ユーザーの承認を得てから**`devtool/flow.sh finish <branch>`を実行する(`develop`へ`--no-ff`マージし、worktreeとブランチを削除する)。マージ後、タスク指示書に「マージ済み(`<merge commit>`)」と記録し、非公開リポジトリにコミットする。
6. **pushしない**: pushはユーザーの指示があったときだけ行う。
7. **まとめる**: 状況を聞かれたら`devtool/tasks.sh`(未マージのTASKと作業ブランチの一覧。`--all`でマージ済みも)で確認して報告する。`devtool/`やドキュメントの自分の変更も`chore/*`ブランチで行い、ユーザーの承認を得てからマージする。
8. **ほかのリポジトリとのやり取り**: ほかのリポジトリ(ProteinDF_bridgeなど)への依頼と返事は、非公開リポジトリの`requests/`にメモとして置く(決まりは`requests/README.md`)。作業を始めるときに、自分宛ての未完了の依頼(`*-to-pytools_*.md`)がないか確認する。対応したら状態を「完了」にし、必要なら返事を書く。

## Claudeがやらないこと

- 公開されるファイル(このリポジトリに追跡されるファイル、コミットメッセージ)に、ローカルのディレクトリのパスを書くこと。ローカルの場所は環境変数や、git管理外の設定ファイル・非公開リポジトリに置く。
- 機能の実装やバグ修正のコードを書くこと(タスク指示書・ドキュメント・`devtool/`の開発用スクリプトは除く)。レビューで見つけた問題は修正依頼として書く。
- 実装担当のworktreeのファイルを編集すること、そこでコミットすること。
- タスク指示書、`TODO.md`、`SPEC.md`をProteinDF_pytools本体にコミット・pushすること(いずれも非公開)。非公開リポジトリをリモートにpushするのは、ユーザーの指示があったときだけ。

## Claudeの作業場所

- Claudeの作業場所は`develop`のworktreeである。実装担当は別のworktreeで作業するので、`develop`でのコミットが実装担当の作業と衝突することはない。
- ただし、コミット前に`git status`と`git branch --show-current`を確認し、`develop`以外にいないこと、自分のもの以外の変更がステージされていないことを確かめる。

## テスト環境についての注意(2026-10-12時点)

- `devtool/check.sh`は`run_test.sh`のdoctest・unittestを、ブランチのworktreeのソースに対して実行する(`PYTHONPATH`で指定)。`PDF_HOME`にProteinDF本体のソースのパスを設定する。未設定の間はtestsがSKIPになる。
- Pythonは`PDFTOOLS_PYTHON`、なければ`.git/check-venv/`のものを使う。`.git/check-venv/`は`devtool/setup-check-venv.sh`で作る。要件: Python 3.12、`numpy`(2.xでよい)、`h5py`、`matplotlib`、`proteindf_bridge` 2026.10.14以上(2026.10.6以前は`Matrix.set_buffer`が壊れておりnumpy 2でも動かない。2026.10.7より前はmsgpackの整数キーを文字列として読むため、`pdfparam.mpac`の読み込みが失敗する)。`proteindf_bridge`はPyPIにないので、`PDFTOOLS_BRIDGE_DIR`にソースのディレクトリかリリースのwheelを指定する。
- **`sitecustomize.py`のあるvenv(他リポジトリの回帰テスト用のvenvなど)を使わない。** 起動時に`proteindf_tools`の不具合を上書きするため、不具合が隠れてテストが通ってしまう。`devtool/setup-check-venv.sh`は、見つけると警告する。
- ローカルの場所は、公開ファイルではなく環境変数や`.git/`内の設定に置く。
