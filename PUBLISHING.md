# Publishing to PyPI

The distribution name is **`para-gsa`**; the import package is **`para_gsa`**.
The version is read from `para_gsa/__init__.py` (`__version__`).

## 1. Local build (sanity check)

```bash
python -m pip install --upgrade build twine
python -m build              # -> dist/para_gsa-<version>.tar.gz + .whl
python -m twine check dist/*
```

The wheel contains only `para_gsa/` (plus `LICENSE`); notebook, tests and the
database-initialization script are not shipped.

## 2. Option A — Trusted Publishing via GitHub Actions (recommended)

No API token is stored anywhere: PyPI trusts the GitHub OIDC identity of this
repository.

1. Register a **pending publisher** on PyPI:
   - <https://pypi.org/manage/account/publishing/>
   - PyPI Project Name: `para-gsa`
   - Owner: `lucas-sys`
   - Repository name: `para_gsa`
   - Workflow name: `publish.yml`
   - Environment name: `pypi`
2. (Optional but recommended) Repeat on TestPyPI, environment `testpypi`:
   - <https://test.pypi.org/manage/account/publishing/>
3. In GitHub, open **Settings → Environments** and create `pypi` (and
   `testpypi`). No secrets are needed.
4. Publish:
   - TestPyPI: **Actions → Publish to PyPI → Run workflow**, target `testpypi`.
   - PyPI: create and push a tag:
     ```bash
     git tag v0.1.0
     git push origin v0.1.0
     ```

## 3. Option B — Manual upload with an API token

1. Create a token at <https://pypi.org/manage/account/token/>.
2. Build the artifacts (section 1), then:
   ```bash
   # TestPyPI first
   python -m twine upload --repository-url https://test.pypi.org/legacy/ dist/*
   # PyPI
   python -m twine upload dist/*
   ```
   Username: `__token__`, password: the token.

## 4. Verify the published package

```bash
# TestPyPI (runtime dependencies come from PyPI)
pip install -i https://test.pypi.org/simple/ --extra-index-url https://pypi.org/simple/ para-gsa

# PyPI
pip install para-gsa
python -c "import para_gsa; print(para_gsa.__version__)"
```

## 5. Release checklist

1. Bump `__version__` in `para_gsa/__init__.py`.
2. Update `README.md` / `README_CN.md` if behaviour changed.
3. Commit, then `git tag vX.Y.Z && git push origin vX.Y.Z`.
4. Check the workflow run and the project page on PyPI.
