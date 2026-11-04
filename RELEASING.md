# Releasing Nanomotif

A release publishes two packages to PyPI from one GitHub release:

- `nanomotif-rs`: the compiled kernels (`rust/`), as prebuilt wheels for Linux x86_64/aarch64, macOS
  (Intel and Apple silicon) and Windows x86_64. No source package is uploaded, so pip never compiles.
- `nanomotif`: the Python package. It depends on `nanomotif-rs` (same version) on those platforms;
  elsewhere it runs its pure-Python code.

`.github/workflows/release.yml` builds the wheels, uploads `nanomotif-rs` (trusted publishing, no
token), then builds and uploads `nanomotif` (with the `PYPI_API_TOKEN` secret). It runs when a GitHub
release is published. Nothing is uploaded unless every wheel builds.

## One-time setup

Needed once, before the first release with `nanomotif-rs` (1.2.1).

1. **GitHub environment.** Repository → Settings → Environments → New environment → name it `pypi`.
   No secrets or rules are needed in it.
2. **PyPI trusted publisher for `nanomotif-rs`.** Log in to pypi.org with an account that should own the
   project → Your account → Publishing → "Add a new pending publisher" → GitHub, and enter:
   - PyPI Project Name: `nanomotif-rs`
   - Owner: `MicrobialDarkMatter`
   - Repository name: `nanomotif`
   - Workflow name: `release.yml`
   - Environment name: `pypi`

   The project is created on the first upload. Afterwards it appears under Your projects; add
   co-maintainers there if wanted.
3. **`PYPI_API_TOKEN` secret** (already used for `nanomotif`): Settings → Secrets and variables →
   Actions. It must be allowed to upload `nanomotif`; it is not used for `nanomotif-rs`.

## Each release

1. **Version.** Set the same version in all four places:
   - `nanomotif/_version.py`
   - `rust/pyproject.toml` (`version`)
   - `rust/Cargo.toml` (`version`)
   - `setup.py` (the `nanomotif-rs==<version>` requirement)

   Add a section to `CHANGELOG.md`.
2. **Tests pass on `main`.** The "Python Build and Test" workflow runs on every push to `main`; check that
   the run for the release commit is green (Actions tab).
3. **Tag and push.**
   ```
   git tag -a v<version> -m "Nanomotif <version>"
   git push origin main v<version>
   ```
4. **Publish the GitHub release.** Releases → Draft a new release → choose tag `v<version>` → title
   `Nanomotif <version>` → paste the version's `CHANGELOG.md` section as the description → Publish release.
5. **Watch the workflow.** Actions → "Release to PyPI". It runs, in order:
   1. `rust-linux`, `rust-macos`, `rust-windows`: build the wheels (about 5-10 minutes)
   2. `publish-rust`: uploads `nanomotif-rs` (wheels only)
   3. `publish-nanomotif`: builds and uploads `nanomotif`

   If a build job fails, nothing has been uploaded: fix, push, and rerun the workflow (Actions → the run →
   "Re-run all jobs", or "Run workflow" on `release.yml` for the tag). If `publish-rust` succeeded and
   `publish-nanomotif` failed, rerun only the failed job; PyPI refuses files that already exist.
6. **Check the result** in a fresh environment, on a machine without a compiler if possible:
   ```
   python -m venv /tmp/nm && /tmp/nm/bin/pip install nanomotif==<version>
   /tmp/nm/bin/python -c "import nanomotif, nanomotif_rs; print('ok')"
   /tmp/nm/bin/nanomotif check_installation
   ```
   pip should download a `nanomotif_rs-<version>-cp39-abi3-<platform>.whl` and build nothing.
7. **Bioconda** (if used): the Bioconda bot opens a version-bump PR for `nanomotif` after the PyPI upload.
   From 1.2.1 the recipe needs `nanomotif-rs` too: either a `nanomotif-rs` recipe built from the GitHub
   release tarball (`rust/`, with `{{ compiler('rust') }}` and `maturin`), or the recipe builds
   `rust/` itself. Without it, the conda package runs the pure-Python code.

## Notes

- A PyPI version can never be uploaded twice. If a release goes wrong after an upload, bump the patch
  version and release again.
- `nanomotif-rs` wheels use the stable Python ABI (abi3, Python >= 3.9): one wheel per platform covers
  all Python versions.
- Without `nanomotif-rs`, Nanomotif gives the same results, but bgzipped pileups are read per bin with
  pysam and plain pileups are loaded whole into memory.
