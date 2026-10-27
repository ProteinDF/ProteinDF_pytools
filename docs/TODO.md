# Documentation TODO

- [ ] Enable GitHub Pages for this repo (Settings → Pages → Source: GitHub Actions)
      so the `docs` workflow can deploy.
- [ ] Translate the API reference and CLI reference tables into Japanese
      (`docs/locale/ja/LC_MESSAGES/api/*.po`, `.../usage.po`'s per-script
      entries were translated, but docstrings inside the modules themselves
      are still English-only) — currently untranslated and falls back to
      English.
- [ ] Install a local TeX distribution to build PDFs
      (`brew install --cask mactex-no-gui`, or `texlive-lang-japanese` on
      Linux for the Japanese PDF) and verify `make latexpdf` /
      `make latexpdf-ja`.
- [ ] Review the auto-generated API reference (`docs/api/*.rst`) for modules
      that should be excluded or reorganized (e.g. abstract base classes
      not meant for direct end-user use).
- [ ] Most `proteindf_tools` modules have little or no docstring content,
      so the autodoc-generated API pages are currently thin. Add
      docstrings to the library's public classes/methods for the API
      reference to be genuinely useful.
- [ ] `proteindf_tools` imports `proteindf_bridge` at runtime but does not
      declare it in `setup.cfg`'s `install_requires`. The docs workflow
      installs it explicitly from git; consider adding it as a proper
      dependency of the package itself.
