"""PubChem REST API wrapper.

Public API
----------
All symbols listed in ``pubchem.pubchem.__all__`` are re-exported at
the package level for convenience::

    from pubchem import cas_to_cid, cid_to_cas, get_from_cids
"""

__version__ = "1.0.0"

from .throttle import safe_request
from .pubchem import (
    PubchemInputError,
    format_cas,
    cas_to_cid,
    cas_to_sid,
    cas_to_pubchem,
    synonyms_to_pubchem,
    inchikey_to_pubchem,
    cas_to_inchi,
    cids_to_cas_and_einecs,
    cids_to_cas_and_einecs_and_dtx,
    single_synonym_to_pubchem,
    SMILES_to_pubchem,
    get_from_cids,
    get_chunks,
    get_mols_from_cas,
    get_cids_from_sids,
    get_mols_from_cids,
    get_mols_from_sids,
    cids_to_mol,
    get_cids_from_smiles,
    cas_to_mols,
    cid_to_cas,
    pubchem_pfas_tree,
    download_cids_by_hnid,
    get_HSDB,
    synonyms_to_smiles,
    synonyms_to_mols,
    dtxsid_to_mols,
    cas_to_smiles,
    dtxsid_to_smiles,
)

__all__ = [
    "__version__",
    "safe_request",
    "PubchemInputError",
    "format_cas",
    "cas_to_cid",
    "cas_to_sid",
    "cas_to_pubchem",
    "synonyms_to_pubchem",
    "inchikey_to_pubchem",
    "cas_to_inchi",
    "cids_to_cas_and_einecs",
    "cids_to_cas_and_einecs_and_dtx",
    "single_synonym_to_pubchem",
    "SMILES_to_pubchem",
    "get_from_cids",
    "get_chunks",
    "get_mols_from_cas",
    "get_cids_from_sids",
    "get_mols_from_cids",
    "get_mols_from_sids",
    "cids_to_mol",
    "get_cids_from_smiles",
    "cas_to_mols",
    "cid_to_cas",
    "pubchem_pfas_tree",
    "download_cids_by_hnid",
    "get_HSDB",
    "synonyms_to_smiles",
    "synonyms_to_mols",
    "cas_to_smiles",
    "dtxsid_to_mols",
    "dtxsid_to_smiles",
]
