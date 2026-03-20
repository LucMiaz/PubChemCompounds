"""
PubChem REST API wrapper utilities.

Provides functions to query the PubChem PUG REST API for compound and
substance information: CID/SID lookups by CAS, InChIKey or SMILES,
batch property fetching, SDF/mol retrieval, and CAS/EINECS extraction.
"""

from __future__ import annotations

from typing import Union
import json
import regex as re
import os
from tqdm import tqdm
import logging
from rdkit import Chem
import pandas as pd
from .throttle import safe_request

logger = logging.getLogger(__name__)

__all__ = [
    "PubchemInputError",
    "format_cas",
    "cas_to_cid",
    "cas_to_sid",
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
]


class PubchemInputError(AttributeError):
    """Raised when PubChem returns an error response for the given input."""


def format_cas(_cas: Union[str, list]) -> str:
    """Ensure a CAS number is formatted as ``XXXXXXXX-YY-Z``.

    Args:
        _cas: CAS number as a string (with or without dashes) or a
            three-element list ``[XXXXXXXX, YY, Z]``.

    Returns:
        Dash-separated CAS string, e.g. ``"7732-18-5"``.
    """
    if (
        isinstance(_cas, list)
        and len(_cas) == 3
        and len(str(_cas[0])) < 8
        and len(str(_cas[1])) == 2
        and len(str(_cas[2])) == 1
    ):
        return f"{_cas[0]}-{_cas[1]}-{_cas[2]}"
    if isinstance(_cas, str) and "-" not in _cas:
        xxxx = _cas[:-3]
        z = _cas[-1]
        yy = _cas[-3:-1]
        if len(xxxx) > 8:
            logger.error("Error with formatting of cas: %s", _cas)
        return f"{xxxx}-{yy}-{z}"
    if isinstance(_cas, str):
        return _cas
    logger.error("Error with formatting of cas : %s", _cas)
    return str(_cas)

def cas_to_cid(cas: Union[str, list]) -> tuple[dict, list]:
    """Convert CAS number(s) to PubChem Compound IDs (CIDs).

    Convenience wrapper around :func:`synonyms_to_pubchem` with
    ``substance=False``.

    Args:
        cas: Single CAS string or a list of CAS strings.

    Returns:
        A tuple ``(mapping, failed)`` where *mapping* is a ``dict``
        keyed by CAS with lists of CIDs as values, and *failed* is a
        list of CAS numbers that could not be resolved.
    """
    return synonyms_to_pubchem(cas, substance=False)


def cas_to_sid(cas: Union[str, list]) -> tuple[dict, list]:
    """Convert CAS number(s) to PubChem Substance IDs (SIDs).

    Convenience wrapper around :func:`synonyms_to_pubchem` with
    ``substance=True``.

    Args:
        cas: Single CAS string or a list of CAS strings.

    Returns:
        A tuple ``(mapping, failed)`` where *mapping* is a ``dict``
        keyed by CAS with lists of SIDs as values, and *failed* is a
        list of CAS numbers that could not be resolved.
    """
    return synonyms_to_pubchem(cas, substance=True)


def inchikey_to_pubchem(inchikey: str) -> Union[list, None]:
    """Convert an InChIKey to a list of PubChem CIDs.

    Args:
        inchikey: Standard 27-character InChIKey string.

    Returns:
        List of CIDs, or ``None`` if not found / on error.
    """
    url = (
        f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/"
        f"{inchikey}/cids/JSON"
    )
    response = safe_request(url)
    try:
        data = json.loads(response)
    except json.JSONDecodeError:
        logger.info("JSON decode error for inchikey %s: %s", inchikey, response)
        return None
    if "Fault" in data:
        logger.info("InChIKey %s got fault: %s", inchikey, data.get("Fault", ""))
        return None
    return data.get("IdentifierList", {}).get("CID", None)


# backward-compatible alias
cas_to_pubchem = None  # defined after synonyms_to_pubchem below


def synonyms_to_pubchem(synonym: Union[str, list], substance: bool) -> tuple[dict, list]:
    """Look up PubChem IDs for one or more synonyms.

    Args:
        synonym: Single synonym string or a list of synonym strings.
        substance: If ``True``, return Substance IDs (SIDs); if
            ``False``, return Compound IDs (CIDs).

    Returns:
        A tuple ``(mapping, failed)`` where *mapping* maps each CAS
        string to its list of IDs, and *failed* contains CAS strings
        that returned no result.
    """
    if synonym is None:
        return {}, []
    if isinstance(synonym, str):
        fsynonym = [synonym]
    else:
        fsynonym = list(synonym)

    return_dict: dict = {}
    failed: list = []
    if len(fsynonym)>1:
        for _synonym in tqdm(fsynonym, total=len(fsynonym)):
            x = single_synonym_to_pubchem(_synonym, substance=substance)
            if x is None:
                failed.append(_synonym)
            else:
                return_dict[_synonym] = x
    else:
        x = single_synonym_to_pubchem(fsynonym[0], substance=substance)
        if x is None:
            failed.append(fsynonym[0])
        else:
            return_dict[fsynonym[0]] = x
    return return_dict, failed


# Backward-compatible alias (was renamed to synonyms_to_pubchem)
cas_to_pubchem = synonyms_to_pubchem


def get_chunks(l: list, n: int) -> list:
    """Divide a list into chunks of at most *n* elements.

    Args:
        l: Input list to split.
        n: Maximum chunk size.

    Returns:
        List of sub-lists, each of length at most *n*.
    """
    try:
        return [l[i : i + n] for i in range(0, len(l), max(1, n))]
    except TypeError as e:
        logger.info("get_chunks error: %s", e)
        return []


def cas_to_inchi(cas: Union[str, list], max_query: int = 100) -> tuple[dict, list]:
    """Convert CAS number(s) to InChI string(s) via PubChem.

    Resolves CAS numbers to CIDs first, then fetches the InChI for
    each CID in batches.

    Args:
        cas: Single CAS string or a list of CAS strings.
        max_query: Maximum number of CIDs per PubChem batch request.

    Returns:
        A tuple ``(cas_inchi, failed)`` where *cas_inchi* maps each
        resolved CAS to its InChI string (or a list of strings if
        multiple CIDs matched), and *failed* lists unresolved CAS
        numbers.
    """
    cas_cids, failed = synonyms_to_pubchem(cas, substance=False)

    cas_inchi: dict = {}
    all_cids = list(set(cid for cids in cas_cids.values() for cid in cids))

    if all_cids:
        cids_chunks = get_chunks(all_cids, max_query)
        cid_to_inchi: dict = {}

        for chunk in cids_chunks:
            properties = get_from_cids(chunk, target="property/InChI")
            if properties:
                for prop in properties.get("PropertyTable", {}).get("Properties", []):
                    cid = prop.get("CID")
                    inchi = prop.get("InChI")
                    if cid and inchi:
                        cid_to_inchi[cid] = inchi

        for cas_num, cids in cas_cids.items():
            inchi_list = [cid_to_inchi[cid] for cid in cids if cid in cid_to_inchi]
            if inchi_list:
                cas_inchi[cas_num] = inchi_list[0] if len(inchi_list) == 1 else inchi_list

    return cas_inchi, failed


def cids_to_cas_and_einecs(cids: list, max_query: int = 100) -> dict:
    """Legacy alias for :func:`cids_to_cas_and_einecs_and_dtx`."""
    return cids_to_cas_and_einecs_and_dtx(cids, max_query)


def cids_to_cas_and_einecs_and_dtx(cids: list, max_query: int = 100) -> dict:
    """Extract CAS, EINECS, and DTXSID identifiers from CID synonyms.

    For each CID, the PubChem synonym list is searched for strings
    matching the CAS (``XXXXXXX-YY-Z``), EINECS (``XXX-XXX-X``), and
    DTXSID (``DTXSIDxxxxxxxx``) patterns.

    Args:
        cids: List of PubChem CIDs.
        max_query: Maximum number of CIDs per batch request.

    Returns:
        Dict mapping each CID to a nested dict with optional keys
        ``"CAS"``, ``"EINECS"``, and ``"DTXSID"``.
    """
    PAT_cas = r"^\d{2,7}[-]\d{2}[-]\d$"
    PAT_einecs = r"^\d{3}[-]\d{3}[-]\d$"
    PAT_dtx = r"(?<=DTXSID)\d{5,10}"
    cids_cas_einecs: dict = {}
    for chunk in tqdm(get_chunks(cids, max_query)):
        synonyms = get_from_cids(chunk, target="synonyms")
        if synonyms is not None:
            for v in synonyms.get("InformationList", {}).get("Information", []):
                for candidate in v.get("Synonym", []):
                    cas_match = re.findall(PAT_cas, candidate)
                    einecs_match = re.findall(PAT_einecs, candidate)
                    dtx_match = re.findall(PAT_dtx, candidate)
                    if cas_match:
                        cids_cas_einecs.setdefault(v["CID"], {})["CAS"] = cas_match[0]
                    if einecs_match:
                        cids_cas_einecs.setdefault(v["CID"], {})["EINECS"] = einecs_match[0]
                    if dtx_match:
                        cids_cas_einecs.setdefault(v["CID"], {})["DTXSID"] = dtx_match[0]
    return cids_cas_einecs

def single_cas_to_pubchem(cas: str, substance: bool = False) -> Union[list, None]:
    """Legacy alias for :func:`single_synonym_to_pubchem`."""
    return single_synonym_to_pubchem(cas, substance)


def single_synonym_to_pubchem(
    synonym: str, substance: bool = False
) -> Union[list, None]:
    """Use the PubChem REST API to convert a synonym to a CID or SID list.

    Args:
        synonym: Any PubChem-recognised synonym (CAS, DTXSID, name, …).
        substance: If ``True``, return Substance IDs (SIDs); if
            ``False`` (default), return Compound IDs (CIDs).

    Returns:
        List of CIDs or SIDs, or ``None`` if the synonym was not found
        or the response was malformed.
    """
    if substance:
        url = (
            f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/substance/name/"
            f"{synonym}/sids/JSON"
        )
    else:
        url = (
            f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/"
            f"{synonym}/cids/JSON"
        )
    response = safe_request(url)
    try:
        data = json.loads(response)
    except json.JSONDecodeError:
        logger.info("JSON decode error for synonym %s: %s", synonym, response)
        return None
    if "Fault" in data:
        logger.info("Synonym %s got fault: %s", synonym, data.get("Fault", ""))
        return None
    if substance:
        return data.get("IdentifierList", {}).get("SID", [])
    return data.get("IdentifierList", {}).get("CID", [])


def SMILES_to_pubchem(smiles: str) -> Union[list, None]:
    """Convert a SMILES string to a list of PubChem CIDs.

    The ``#`` character (triple bond) is URL-encoded automatically.

    Args:
        smiles: SMILES string of the compound.

    Returns:
        List of CIDs, or ``None`` if not found or on error.
    """
    encoded = smiles.replace("#", "%23")
    url = (
        f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/"
        f"{encoded}/cids/JSON"
    )
    response = safe_request(url)
    try:
        data = json.loads(response)
    except json.JSONDecodeError:
        logger.info("JSON decode error for SMILES %s: %s", smiles, response)
        return None
    if "Fault" in data:
        logger.info("SMILES %s got fault: %s", smiles, data.get("Fault", ""))
        return None
    return data.get("IdentifierList", {}).get("CID", [])


def get_from_cids(
    cids: Union[list, int],
    target: Union[str, None] = None,
    format: str = "JSON",  # noqa: A002
) -> Union[dict, None]:
    """Fetch data from PubChem for the given CIDs.

    Args:
        cids: Single CID (int) or list of CIDs.
        target: PubChem endpoint suffix, e.g. ``"synonyms"`` or
            ``"property/InChI"``.  When ``None``, the compound JSON
            record itself is retrieved.
        format: Response format (default: ``"JSON"``).

    Returns:
        Parsed JSON response as a dict, or ``None`` on error.
    """
    if cids is not None and not isinstance(cids, list):
        cids = [cids]
    fcid = ",".join(str(x) for x in cids)
    if target is not None:
        url = (
            f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/"
            f"{fcid}/{target}/{format}"
        )
    else:
        url = (
            f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/"
            f"{fcid}/{format}"
        )
    response = safe_request(url)
    try:
        return json.loads(response)
    except Exception as e:  # pylint: disable=broad-except
        logger.error("Could not fetch %s: %s", url, e)
        return None


def get_mols_from_cas(
    cas: Union[list, str], substance: bool = False
) -> Chem.SDMolSupplier:
    """Fetch RDKit molecules for CAS number(s) via PubChem.

    Args:
        cas: Single CAS string or list of CAS strings.
        substance: If ``True``, resolve via SIDs instead of CIDs.

    Returns:
        An :class:`rdkit.Chem.SDMolSupplier` for the matched structures.
    """
    cids, _ = synonyms_to_pubchem(cas, substance=substance)
    return get_mols_from_cids(list(cids.values()))


def get_cids_from_sids(
    sids: Union[list, int], max_query: int = 50
) -> dict:
    """Map PubChem Substance IDs (SIDs) to Compound IDs (CIDs).

    Args:
        sids: Single SID or list of SIDs.
        max_query: Maximum number of SIDs per batch request.

    Returns:
        Dict mapping each SID to its CID (or ``None`` if unmapped).
    """
    if sids is not None and not isinstance(sids, list):
        sids = [sids]
    ret: dict = {}
    for csids in get_chunks(sids, max_query):
        fsid = ",".join(str(x) for x in csids)
        url = (
            f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/substance/sid/"
            f"{fsid}/cids/JSON"
        )
        raw = safe_request(url)
        try:
            data = json.loads(raw)
        except Exception as e:  # pylint: disable=broad-except
            logger.error("get_cids_from_sids JSON error: %s", e)
            raise
        for info in data.get("InformationList", {}).get("Information", []):
            ret[info["SID"]] = info.get("CID", None)
    return ret


def get_mols_from_cids(
    cids: Union[list, int],
    index: int = 0,
    max_query: int = 500,
    filename: Union[str, None] = None,
) -> tuple:
    """Fetch SDF data from PubChem and return an RDKit mol supplier.

    Downloaded SDF chunks are appended to *filename* and the file is
    kept on disk so the caller can reuse or delete it.

    Args:
        cids: Single CID or list of CIDs.
        index: Suffix for the temporary filename when *filename* is
            not provided.
        max_query: Maximum number of CIDs per batch request.
        filename: Path for the output SDF file.  A temporary name is
            used when ``None``.

    Returns:
        A tuple ``(supplier, filename)`` where *supplier* is an
        :class:`rdkit.Chem.SDMolSupplier` (or an empty list on error)
        and *filename* is the path of the written SDF file.

    Raises:
        PubchemInputError: If the server returns an XML error page.
    """
    if cids is not None and not isinstance(cids, list):
        cids = [cids]
    if filename is None:
        filename = f"sdffrompubchem_temp_{index}.sdf"
    for ccids in get_chunks(cids, max_query):
        fcid = ",".join(str(x) for x in ccids)
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{fcid}/SDF"
        data = safe_request(url)
        if data[:5] == b"<?xml":
            logger.error("get_mols_from_cids: server returned an XML error page")
            raise PubchemInputError(
                "get_mols_from_cids: received an error from server"
            )
        with open(filename, "ab") as f:
            f.write(data)
    try:
        suppl = Chem.SDMolSupplier(filename, sanitize=True, removeHs=False)
        return suppl, filename
    except Exception as e:  # pylint: disable=broad-except
        logger.info("get_mols_from_cids error for index %s: %s", index, e)
        return [], ""


def get_mols_from_sids(
    sids: Union[list, int],
    index: int = 0,
    max_query: int = 200,
    filename: Union[str, None] = None,
) -> tuple:
    """Fetch SDF data from PubChem for substance IDs (SIDs).

    Analogous to :func:`get_mols_from_cids` but queries the substance
    endpoint.

    Args:
        sids: Single SID or list of SIDs.
        index: Suffix for the temporary filename when *filename* is
            not provided.
        max_query: Maximum number of SIDs per batch request.
        filename: Path for the output SDF file.

    Returns:
        A tuple ``(supplier, filename)`` — see :func:`get_mols_from_cids`.

    Raises:
        PubchemInputError: If the server returns an XML error page.
    """
    if sids is not None and not isinstance(sids, list):
        sids = [sids]
    if filename is None:
        filename = f"sdffrompubchem_temp_{index}.sdf"
    for ssids in get_chunks(sids, max_query):
        fsid = ",".join(str(x) for x in ssids)
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/substance/sid/{fsid}/SDF"
        data = safe_request(url)
        if data[:5] == b"<?xml":
            logger.error("get_mols_from_sids: server returned an XML error page")
            raise PubchemInputError(
                "get_mols_from_sids: received an error from server"
            )
        with open(filename, "ab") as f:
            f.write(data)
    try:
        suppl = Chem.SDMolSupplier(filename, sanitize=True, removeHs=False)
        return suppl, filename
    except Exception as e:  # pylint: disable=broad-except
        logger.info("get_mols_from_sids error for index %s: %s", index, e)
        return [], ""

def cids_to_mol(
    cids: Union[list, int],
    filename: Union[str, None] = None,
    max_query: int = 500,
):
    """Yield RDKit molecules for the given CIDs.

    Args:
        cids: Single CID or list of CIDs.
        filename: Optional path for the intermediate SDF file.
        max_query: Maximum number of CIDs per batch request.

    Yields:
        :class:`rdkit.Chem.Mol` instances (may include ``None`` for
        unparseable entries).
    """
    suppl, _ = get_mols_from_cids(cids, filename=filename, max_query=max_query)
    for mol in suppl:
        yield mol


def get_cids_from_smiles(smiles: str) -> Union[list, None]:
    """Convert a SMILES string to a list of PubChem CIDs.

    Args:
        smiles: SMILES string of the compound.

    Returns:
        List of CIDs, or ``None`` if not found.
    """
    url = (
        f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/"
        f"{smiles}/cids/JSON"
    )
    response = safe_request(url)
    data = json.loads(response)
    return data.get("IdentifierList", {}).get("CID", None)


def synonyms_to_mols(
    cas: Union[list, str],
    save: Union[str, None] = None,
    max_cids: Union[int, None] = None,
) -> dict:
    """Fetch RDKit molecules for synonym(s) or other identifiers.

    Resolves each synonym to CIDs via PubChem, downloads the
    corresponding SDF, sets the ``"CAS"`` property on every molecule,
    and returns a dict keyed by synonym.

    Args:
        synonyms: Single synonym string (CAS, DTXSID, …) or list of
            synonym strings.
        save: If provided, the synonym→CID mapping is saved to
            ``<save>_cas_cids.json`` before fetching structures.
        max_cids: Maximum number of CIDs to fetch per synonym.  When
            set to ``1`` only the canonical/preferred CID
            returned by PubChem is used, which prevents hundreds of
            stereoisomers or salt forms from being returned for a
            single DTXSID or CAS number.  Set to ``None`` (default) to fetch
            all matching CIDs.

    Returns:
        Dict mapping each synonym to a list of
        :class:`rdkit.Chem.Mol` objects.
    """
    synonyms_cids, _ = synonyms_to_pubchem(cas, substance=False)
    if save is not None:
        with open(save + "_cas_cids.json", "w", encoding="utf-8") as f:
            json.dump(synonyms_cids, f)
    mols: dict = {}
    files: list = []
    for i, (cas_num, cids) in enumerate(synonyms_cids.items()):
        if max_cids is not None:
            cids = cids[:max_cids]
        suppl, filename = get_mols_from_cids(cids, index=i)
        for m in suppl:
            if m is not None:
                m.SetProp("CAS", cas_num)
                mols.setdefault(cas_num, []).append(m)
        files.append(filename)
    for filename in files:
        try:
            os.remove(filename)
        except (FileNotFoundError, PermissionError):
            pass
    return mols

def cas_to_mols(cas: Union[list, str], max_cids: Union[int, None] = None) -> dict:
    """Fetch RDKit molecules for CAS number(s) via PubChem.
_
    Args:
        cas: Single CAS string or list of CAS strings.
        max_cids: Maximum number of CIDs to fetch per CAS.  When set
            to ``1``, only the canonical/preferred CID returned by
            PubChem is used, which prevents hundreds of stereoisomers
            or salt forms from being returned for a single CAS number.
            Set to ``None`` (default) to fetch all matching CIDs.
    """
    return synonyms_to_mols(cas, max_cids=max_cids)

def cid_to_cas(
    cid: Union[str, int, list], get_einecs: bool = False
) -> Union[dict, list, None]:
    """Convert PubChem CID(s) to CAS numbers (and optionally EINECS).

    Args:
        cid: Single CID (int or str) or a list of CIDs.
        get_einecs: If ``True``, also return EINECS numbers.

    Returns:
        When *get_einecs* is ``False``:

        * Single CID → list of CAS strings for that CID.
        * Multiple CIDs → ``dict`` mapping each CID to its CAS list.

        When *get_einecs* is ``True``:

        * Dict mapping each CID to a ``(cas_list, einecs_list)`` tuple.
    """
    if isinstance(cid, (str, int)):
        multiple = False
        cid = [cid]
    elif isinstance(cid, list):
        multiple = True
    else:
        logger.info("cid_to_cas: %s has an invalid format", cid)
        return None

    data = get_from_cids(cid, "synonyms")
    if data is None:
        return None

    cas: dict = {}
    einecs: dict = {}
    PAT = r"\d{1,10}[-]\d{2}[-]\d"
    EINECS_PAT = r"(EINECS\s)?(\d{3}[-]\d{3}[-]\d)"

    for info in data.get("InformationList", {}).get("Information", []):
        icid = info.get("CID")
        if icid is None:
            continue
        for synonym in info.get("Synonym", []):
            if re.fullmatch(PAT, synonym):
                cas.setdefault(icid, []).append(synonym)
            if get_einecs:
                match = re.fullmatch(EINECS_PAT, synonym)
                if match:
                    einecs.setdefault(icid, []).append(match.group(2))

    if not get_einecs:
        return cas if multiple else cas.get(cid[0], [None])

    combined = {
        c: (cas.get(c, [None]), einecs.get(c, [None]))
        for c in set(cas) | set(einecs)
    }
    return combined if multiple else combined.get(cid[0], None)


def pubchem_pfas_tree(hnid: int = 5517102, max_tries: int = 2) -> Union[list, None]:
    """Fetch CIDs belonging to a PubChem classification tree node.

    Args:
        hnid: Hierarchy node ID.  Defaults to the OECD PFAS list node
            (5517102).
        max_tries: Unused; kept for backwards compatibility.

    Returns:
        List of CIDs in the hierarchy, or ``None`` on error.
    """
    url = (
        f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/classification/hnid/"
        f"{hnid}/cids/json"
    )
    response = safe_request(url)
    try:
        data = json.loads(response)
    except Exception:  # pylint: disable=broad-except
        return None
    id_list = data.get("IdentifierList")
    if id_list is None:
        return None
    return id_list.get("CID")


def download_cids_by_hnid(
    hnid: int, folder_path: str, force: bool = False
) -> Union[list, None]:
    """Download and cache the CID list for a PubChem hierarchy node.

    The result is stored in ``<folder_path>/hnid_XXXXXXX.json``.  On
    subsequent calls the cached file is read unless *force* is
    ``True``.

    Args:
        hnid: Hierarchy node ID.
        folder_path: Directory where the JSON cache files are stored.
        force: If ``True``, always re-download even when a cache file
            already exists.

    Returns:
        List of CIDs, or ``None`` if the download failed.
    """
    filename = f"{folder_path}/hnid_{str(hnid).zfill(7)}.json"
    if force or not os.path.exists(filename):
        cids_all = pubchem_pfas_tree(hnid=hnid)
        with open(filename, "w", encoding="utf-8") as f:
            json.dump(cids_all, f)
        return cids_all
    with open(filename, "r", encoding="utf-8") as f:
        return json.load(f)


def get_HSDB(cids: Union[list, int]) -> Union[dict, int, None]:
    """Return HSDB identifier(s) for the given CID(s).

    Args:
        cids: Single CID or list of CIDs.

    Returns:
        When a single CID is given, returns the HSDB number (int) or
        ``None``.  For a list of CIDs, returns a dict mapping each CID
        to its HSDB number (or ``None``).

    Raises:
        ValueError: If *cids* is not an int or a list.
    """
    PAT = r"(?<=HSDB\s)(\d*)"

    def _find_hsdb(synonyms: list) -> Union[int, None]:
        for synonym in synonyms:
            for match in re.finditer(PAT, synonym):
                return int(match.group())
        return None

    if isinstance(cids, int):
        cids = [cids]
    elif not isinstance(cids, list):
        raise ValueError("cids must be an int or a list of ints")
    info = (
        get_from_cids(cids, target="synonyms")
        .get("InformationList", {})
        .get("Information", {})
    )
    data = {d["CID"]: _find_hsdb(d["Synonym"]) for d in info}
    if len(cids) == 1:
        return data[cids[0]]
    return data


# ---------------------------------------------------------------------------
# SMILES conversion helpers
# ---------------------------------------------------------------------------


def synonyms_to_smiles(
    synonym_list: Union[list, str],
    unique: bool = True,
    join_smiles: bool = True,
) -> tuple[dict, list]:
    """Convert a list of synonyms (CAS, DTXSID, …) to SMILES strings.

    For each synonym the corresponding molecules are fetched from
    PubChem.  When *unique* is ``True`` duplicates are removed by
    comparing InChIKeys.  Multiple structures per synonym are joined
    with a ``.`` when *join_smiles* is ``True``, or returned as a list
    otherwise.

    Only the canonical (first) CID returned by PubChem is used per
    synonym by default; pass ``max_cids=None`` to :func:`cas_to_mols`
    directly if you need all forms.

    Args:
        synonym_list: Synonym or list of synonyms to resolve (CAS
            numbers, DTXSIDs, compound names, …).  A bare string is
            treated as a single-element list.
        unique: Deduplicate structures by InChIKey.
        join_smiles: If ``True``, join multiple SMILES with ``"."``
            into a single string; otherwise return a list.

    Returns:
        A tuple ``(processed, failed)`` where *processed* maps each
        synonym to its SMILES (str or list), and *failed* contains
        synonyms that could not be resolved.
    """
    if isinstance(synonym_list, str):
        synonym_list = [synonym_list]

    processed: dict = {}
    failed: list = []

    for synonym in synonym_list:
        if not synonym:
            continue
        try:
            mols_dict = cas_to_mols(synonym)
            if synonym not in mols_dict:
                failed.append(synonym)
                processed[synonym] = None
                continue
            mol_list = mols_dict[synonym]
            smiles_list: list = []
            if unique:
                seen_inchikeys: set = set()
                for mol in mol_list:
                    if mol is None:
                        continue
                    try:
                        inchi_key = Chem.MolToInchiKey(mol)
                        if inchi_key not in seen_inchikeys:
                            smiles_list.append(Chem.MolToSmiles(mol))
                            seen_inchikeys.add(inchi_key)
                    except Exception as e:  # pylint: disable=broad-except
                        logger.warning(
                            "Could not process molecule for %s: %s", synonym, e
                        )
            else:
                smiles_list = [
                    Chem.MolToSmiles(m) for m in mol_list if m is not None
                ]

            if smiles_list:
                processed[synonym] = (
                    ".".join(smiles_list) if join_smiles else smiles_list
                )
            else:
                failed.append(synonym)
                processed[synonym] = None
        except Exception as e:  # pylint: disable=broad-except
            logger.error("Error processing synonym %s: %s", synonym, e)
            failed.append(synonym)
            processed[synonym] = None

    return processed, failed


def cas_to_smiles(
    cas_list: list,
    unique: bool = True,
    join_smiles: bool = True,
) -> tuple[dict, list]:
    """Convert CAS numbers to SMILES strings.

    Convenience wrapper around :func:`synonyms_to_smiles`.

    Args:
        cas_list: List of CAS numbers.
        unique: Deduplicate structures by InChIKey.
        join_smiles: Join multiple SMILES per CAS with ``"."``.

    Returns:
        ``(processed, failed)`` — see :func:`synonyms_to_smiles`.
    """
    return synonyms_to_smiles(cas_list, unique=unique, join_smiles=join_smiles)


def dtxsid_to_smiles(
    dtxsid_list: list,
    unique: bool = True,
    join_smiles: bool = True,
) -> tuple[dict, list]:
    """Convert DTXSID identifiers (or any PubChem synonyms) to SMILES.

    Convenience wrapper around :func:`synonyms_to_smiles`.

    Args:
        dtxsid_list: List of DTXSID strings.
        unique: Deduplicate structures by InChIKey.
        join_smiles: Join multiple SMILES per DTXSID with ``"."``.

    Returns:
        ``(processed, failed)`` — see :func:`synonyms_to_smiles`.
    """
    return synonyms_to_smiles(dtxsid_list, unique=unique, join_smiles=join_smiles)

def dtxsid_to_mols(dtxsid_list: list) -> tuple[dict, list]:
    """Convert DTXSID identifiers to RDKit molecules.

    Args:
        dtxsid_list: List of DTXSID strings.

    Returns:
        A tuple ``(processed, failed)`` where *processed* maps each
        DTXSID to a list of :class:`rdkit.Chem.Mol` objects, and
        *failed* contains DTXSIDs that could not be resolved.
    """
    processed: dict = {}
    failed: list = []

    for dtxsid in dtxsid_list:
        if not dtxsid:
            continue
        try:
            mols_dict = cas_to_mols(dtxsid)
            if dtxsid not in mols_dict:
                failed.append(dtxsid)
                processed[dtxsid] = None
                continue
            processed[dtxsid] = mols_dict[dtxsid]
        except Exception as e:  # pylint: disable=broad-except
            logger.error("Error processing DTXSID %s: %s", dtxsid, e)
            failed.append(dtxsid)
            processed[dtxsid] = None

    return processed, failed