Changelog
=========

1.1.0 (2026-04-03)
-------------------

* Renamed importable package from ``pubchem`` to ``pubchem_compounds``
  (PyPI distribution name: ``pubchem-compounds``).
* Changed licence from LGPLv2+ to **CC BY-NC 4.0**.
* Added ``LICENSE`` file to repository root.
* Updated ``pyproject.toml``: new package name, ``license = { file = "LICENSE" }``,
  removed LGPLv2+ classifier, added ``[project.urls]`` and
  ``[tool.setuptools.packages.find]``, added Python 3.12 classifier.
* Removed legacy ``setup.py`` (fully superseded by ``pyproject.toml``).
* Fixed ``synonyms_to_smiles``: added missing ``max_cids`` parameter to
  the function signature (was already used in the body — silent ``NameError`` at runtime).
* Fixed ``throttle.py``: added ``@wraps`` decorator to the inner wrapper
  closure; rewrote all docstrings to Google style for consistency.
* Removed unused ``import pandas as pd`` from ``pubchem_compounds.py``.
* Updated Sphinx docs: ``currentmodule`` updated, fixed ``synonym_to_smiles``
  (singular) typo → ``synonyms_to_smiles`` (plural) in ``api.rst``.
* Added ``docs/source/examples.rst`` with practical usage examples.

1.0.0 (2026-03-19)
-------------------

* Complete rewrite with consistent Google-style docstrings on all public functions.
* Added ``__all__`` to both ``pubchem.pubchem`` and the package ``__init__``;
  all public symbols are now importable directly from ``pubchem``.
* Fixed bytes comparison bug in :func:`~pubchem_compounds.get_mols_from_cids` and
  :func:`~pubchem_compounds.get_mols_from_sids` (``data[:5] == b"<?xml"`` instead of
  comparing bytes to a plain string).
* Removed unused imports (``pandas``, ``BytesIO``, bare ``logging`` calls
  replaced by module-level :mod:`logging` logger).
* :func:`~pubchem_compounds.synonyms_to_smiles` fixed to iterate over the correct
  parameter name; added structured error handling and logging.
* Removed all commented-out dead code.
* Added Sphinx documentation (``docs/``), ReadTheDocs configuration, and
  GitHub Actions CI (pylint + pytest).
* Bumped ``requires-python`` to ``>=3.9``.

0.0.4
-----

* Initial public release on GitLab.
