# Contributing to moldrug

## Running the checks

```bash
pip install -e .[dev]
pytest tests                 # the suite downloads AutoDock-Vina on the first run
ruff check .                 # style and docstrings, configured in pyproject.toml
```

`pyproject.toml` holds the configuration of `ruff`, and the same convention under `[tool.pydocstyle]` for anyone running `pydocstyle` or `prospector` directly. The package and the tests are expected to be clean. `Contrib` and `streamlit` are only checked for the `E`, `F` and `W` rules, since they are scripts rather than part of the published API.

Docstrings follow the numpy convention, with the summary on the first line of the docstring.

## Comments that silence a checker

A few lines carry a marker that turns off a warning. They are deliberate, and the reason is next to them, so please leave them in place:

- `# nosec B403` and `# nosec B301` in `utils.py` and `compat.py`. A checkpoint holds the cost function the user wrote, so it has to be unpickled, and unpickling runs whatever the file names. The docstrings of `loosen` and `decompress_pickle` say so.
- `# hadolint ignore=DL3008` in the `Dockerfile`, where the apt packages are deliberately unpinned.

Pull requests are also analysed by [Codacy](https://app.codacy.com/gh/ale94mleon/moldrug). Its settings are described in `.github/MAINTAINING.md`.
