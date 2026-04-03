Installation
============

Using pip
---------

.. code-block:: bash

   pip install pubchem-compounds

RDKit (optional, for mol/SMILES functions)
------------------------------------------

Functions that return RDKit molecules or SMILES strings require RDKit:

.. code-block:: bash

   pip install rdkit
   # or via conda:
   conda install -c conda-forge rdkit

From source
-----------

.. code-block:: bash

   git clone https://gitlab.com/lucmiaz/pubchem.git
   cd pubchem
   pip install -e .

Optional dependencies
---------------------

.. code-block:: bash

   pip install -e ".[dev]"   # pytest, pylint
   pip install -e ".[docs]"  # sphinx, sphinx-rtd-theme

Dependencies
------------

- `requests <https://requests.readthedocs.io>`_
- `numpy <https://numpy.org>`_
- `regex <https://pypi.org/project/regex/>`_
- `tqdm <https://tqdm.github.io>`_
- `rdkit <https://www.rdkit.org>`_ (for SDF/mol functionality)
