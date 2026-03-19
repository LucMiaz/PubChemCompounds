API Reference
=============

.. currentmodule:: pubchem

Core lookup functions
---------------------

.. autofunction:: cas_to_cid
.. autofunction:: cas_to_sid
.. autofunction:: cas_to_pubchem
.. autofunction:: single_synonym_to_pubchem
.. autofunction:: inchikey_to_pubchem
.. autofunction:: SMILES_to_pubchem

Batch CID utilities
-------------------

.. autofunction:: get_from_cids
.. autofunction:: get_chunks
.. autofunction:: cids_to_cas_and_einecs_and_dtx
.. autofunction:: cids_to_cas_and_einecs
.. autofunction:: cid_to_cas

InChI / SMILES conversion
--------------------------

.. autofunction:: cas_to_inchi
.. autofunction:: cas_to_smiles
.. autofunction:: dtxsid_to_smiles
.. autofunction:: synonym_to_smiles

Mol / SDF retrieval
-------------------

.. autofunction:: get_mols_from_cids
.. autofunction:: get_mols_from_sids
.. autofunction:: get_mols_from_cas
.. autofunction:: cas_to_mols
.. autofunction:: cids_to_mol
.. autofunction:: get_cids_from_smiles

SID utilities
-------------

.. autofunction:: get_cids_from_sids

PFAS / classification tree
--------------------------

.. autofunction:: pubchem_pfas_tree
.. autofunction:: download_cids_by_hnid
.. autofunction:: get_HSDB

CAS formatting
--------------

.. autofunction:: format_cas

Exceptions
----------

.. autoexception:: PubchemInputError

Throttle / request utilities
-----------------------------

.. autofunction:: safe_request
.. autoclass:: pubchem.throttle.metered_request_decorator
   :members:
