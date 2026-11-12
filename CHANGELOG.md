# Changelog

All notable changes to pubchem-compounds are listed here, newest first.
Versions follow [semantic versioning](https://semver.org/); dates are those of
the version bump in `pyproject.toml`.

## 1.1.3 (2026-07-04)

- Added `NoStructureError` (a subclass of `PubchemInputError`), raised by `mol_from_comptox` when CompTox has no structure on file for a DTXSID or returns content RDKit cannot parse.
- `mols_from_comptox` skips DTXSIDs without a structure instead of aborting the batch, prints a one-line summary of the skipped ones, and returns `(mols, missing)` with `return_missing=True`.
- The Windows file-lock retry in `mol_from_comptox` now only retries on the "used by another process" error.

## 1.1.2 (2026-07-03)

- `mol_from_comptox` retries reading the downloaded file when Windows antivirus scanning briefly locks it.
- `mol_from_comptox` and `mols_from_comptox` exported at package level.

## 1.1.1 (2026-07-03)

- Added `mol_from_comptox` and `mols_from_comptox`: fetch RDKit molecules for DTXSIDs directly from the EPA CompTox Dashboard.

## 1.1.0 (2026-04-03)

- Renamed importable package from `pubchem` to `pubchem_compounds` (PyPI distribution name: `pubchem-compounds`).
- Changed licence from LGPLv2+ to **CC BY-NC 4.0**.
- Added `LICENSE` file to repository root.
- Updated `pyproject.toml`: new package name, `license = { file = "LICENSE" }`, removed LGPLv2+ classifier, added `[project.urls]` and `[tool.setuptools.packages.find]`, added Python 3.12 classifier.
- Removed legacy `setup.py` (fully superseded by `pyproject.toml`).
- Fixed `synonyms_to_smiles`: added missing `max_cids` parameter to the function signature (was already used in the body — silent `NameError` at runtime).
- Fixed `throttle.py`: added `@wraps` decorator to the inner wrapper closure; rewrote all docstrings to Google style for consistency.
- Removed unused `import pandas as pd` from `pubchem_compounds.py`.
- Updated Sphinx docs: `currentmodule` updated, fixed `synonym_to_smiles` (singular) typo → `synonyms_to_smiles` (plural) in `api.rst`.
- Added `docs/source/examples.rst` with practical usage examples.

## 1.0.0 (2026-03-19)

- Complete rewrite with consistent Google-style docstrings on all public functions.
- Added `__all__` to both `pubchem.pubchem` and the package `__init__`; all public symbols are now importable directly from `pubchem`.
- Fixed bytes comparison bug in `get_mols_from_cids` and `get_mols_from_sids` (`data[:5] == b"<?xml"` instead of comparing bytes to a plain string).
- Removed unused imports (`pandas`, `BytesIO`, bare `logging` calls replaced by module-level `logging` logger).
- `synonyms_to_smiles` fixed to iterate over the correct parameter name; added structured error handling and logging.
- Removed all commented-out dead code.
- Added Sphinx documentation (`docs/`), ReadTheDocs configuration, and GitHub Actions CI (pylint + pytest).
- Bumped `requires-python` to `>=3.9`.
- Also since 0.0.4 (2024–2025, not released separately): chunked queries and missing-data handling in `get_HSDB`; `cas_to_pubchem` rewritten on top of `single_synonym_to_pubchem`, with a tqdm progress bar; `cas_to_inchi` and `cas_to_smiles`; `cas_to_mols` adapted to the current PUG REST API; Python 3.13 support.

## 0.0.4 (2024-08-29)

- `inchikey_to_pubchem` and `SMILES_to_pubchem` lookups; DTXSIDs included in the synonyms fetched for a compound.
- `get_HSDB` to fetch the Hazardous Substances Data Bank records of CIDs.
- Progress bar for CAS lookups; updated synonym fetch.

## 0.0.3 (2024-05-23)

- Molecules from SIDs (`get_mols_from_sids`).
- CAS numbers saved in the SDF written by `cas_to_mols`, which can also save its progress midway.
- CIDs from the PubChem PFAS classification tree (`pubchem_pfas_tree`).
- Failed fetches no longer raise; they are reported instead.

## 0.0.2 (2024-02-21)

- First updates after the initial version.

## 0.0.1 (2024-02-15)

- Initial version, on GitLab (`lucmiaz/pubchem`).
