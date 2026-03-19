Changelog
=========

1.0.0 (2026-03-19)
-------------------

* Complete rewrite with consistent Google-style docstrings on all public functions.
* Added ``__all__`` to both ``pubchem.pubchem`` and the package ``__init__``;
  all public symbols are now importable directly from ``pubchem``.
* Fixed bytes comparison bug in :func:`~pubchem.get_mols_from_cids` and
  :func:`~pubchem.get_mols_from_sids` (``data[:5] == b"<?xml"`` instead of
  comparing bytes to a plain string).
* Removed unused imports (``pandas``, ``BytesIO``, bare ``logging`` calls
  replaced by module-level :mod:`logging` logger).
* :func:`~pubchem.synonym_to_smiles` fixed to iterate over the correct
  parameter name; added structured error handling and logging.
* Removed all commented-out dead code.
* Added Sphinx documentation (``docs/``), ReadTheDocs configuration, and
  GitHub Actions CI (pylint + pytest).
* Bumped ``requires-python`` to ``>=3.9``.

0.0.4
-----

* Initial public release on GitLab.
