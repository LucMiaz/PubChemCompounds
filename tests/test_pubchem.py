"""Tests for the pubchem package.

All tests that would make real HTTP calls mock ``safe_request`` so the
test suite runs without network access.
"""

from __future__ import annotations

import json
from unittest.mock import MagicMock, patch

import pytest

import pubchem
from pubchem import (
    PubchemInputError,
    cas_to_cid,
    cas_to_inchi,
    cas_to_pubchem,
    cas_to_sid,
    cids_to_cas_and_einecs_and_dtx,
    format_cas,
    get_chunks,
    get_cids_from_sids,
    get_from_cids,
    inchikey_to_pubchem,
    single_synonym_to_pubchem,
    SMILES_to_pubchem,
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _json_bytes(payload: dict) -> bytes:
    """Encode *payload* as UTF-8 JSON bytes (mimicking safe_request output)."""
    return json.dumps(payload).encode()


# ---------------------------------------------------------------------------
# Package-level checks
# ---------------------------------------------------------------------------


def test_public_api_exported():
    """All key symbols must be accessible directly on the package."""
    for name in [
        "cas_to_cid",
        "cas_to_sid",
        "cas_to_pubchem",
        "inchikey_to_pubchem",
        "cas_to_inchi",
        "cids_to_cas_and_einecs_and_dtx",
        "single_synonym_to_pubchem",
        "SMILES_to_pubchem",
        "get_from_cids",
        "get_chunks",
        "cid_to_cas",
        "PubchemInputError",
    ]:
        assert hasattr(pubchem, name), f"pubchem.{name} is missing"


# ---------------------------------------------------------------------------
# format_cas
# ---------------------------------------------------------------------------


class TestFormatCas:
    def test_already_formatted(self):
        assert format_cas("7732-18-5") == "7732-18-5"

    def test_no_dashes(self):
        assert format_cas("7732185") == "7732-18-5"

    def test_list_input(self):
        assert format_cas([7732, "18", "5"]) == "7732-18-5"

    def test_long_cas(self):
        # 10 characters before the check digits
        result = format_cas("1234567890-12-3")
        assert result == "1234567890-12-3"


# ---------------------------------------------------------------------------
# get_chunks
# ---------------------------------------------------------------------------


class TestGetChunks:
    def test_even_split(self):
        chunks = get_chunks([1, 2, 3, 4], 2)
        assert chunks == [[1, 2], [3, 4]]

    def test_remainder(self):
        chunks = get_chunks([1, 2, 3, 4, 5], 2)
        assert chunks == [[1, 2], [3, 4], [5]]

    def test_single_chunk(self):
        assert get_chunks([1, 2, 3], 10) == [[1, 2, 3]]

    def test_empty_list(self):
        assert get_chunks([], 5) == []

    def test_chunk_size_one(self):
        assert get_chunks([1, 2, 3], 1) == [[1], [2], [3]]

    def test_invalid_input_returns_empty(self):
        # TypeError from None should be caught
        assert get_chunks(None, 5) == []  # type: ignore[arg-type]


# ---------------------------------------------------------------------------
# single_synonym_to_pubchem
# ---------------------------------------------------------------------------


class TestSingleSynonymToPubchem:
    @patch("pubchem.pubchem.safe_request")
    def test_returns_cids(self, mock_req):
        mock_req.return_value = _json_bytes(
            {"IdentifierList": {"CID": [962]}}
        )
        result = single_synonym_to_pubchem("7732-18-5", substance=False)
        assert result == [962]

    @patch("pubchem.pubchem.safe_request")
    def test_returns_sids(self, mock_req):
        mock_req.return_value = _json_bytes(
            {"IdentifierList": {"SID": [12345]}}
        )
        result = single_synonym_to_pubchem("7732-18-5", substance=True)
        assert result == [12345]

    @patch("pubchem.pubchem.safe_request")
    def test_fault_returns_none(self, mock_req):
        mock_req.return_value = _json_bytes(
            {"Fault": {"Code": "PUGREST.NotFound", "Message": "No CID found"}}
        )
        assert single_synonym_to_pubchem("INVALIDCAS") is None

    @patch("pubchem.pubchem.safe_request")
    def test_invalid_json_returns_none(self, mock_req):
        mock_req.return_value = b"not json"
        assert single_synonym_to_pubchem("7732-18-5") is None


# ---------------------------------------------------------------------------
# cas_to_pubchem / cas_to_cid / cas_to_sid
# ---------------------------------------------------------------------------


class TestCasToPubchem:
    @patch("pubchem.pubchem.safe_request")
    def test_single_cas_cid(self, mock_req):
        mock_req.return_value = _json_bytes(
            {"IdentifierList": {"CID": [962]}}
        )
        mapping, failed = cas_to_pubchem("7732-18-5", substance=False)
        assert mapping == {"7732-18-5": [962]}
        assert failed == []

    @patch("pubchem.pubchem.safe_request")
    def test_multiple_cas(self, mock_req):
        mock_req.return_value = _json_bytes(
            {"IdentifierList": {"CID": [962]}}
        )
        mapping, failed = cas_to_pubchem(["7732-18-5", "74-82-8"], substance=False)
        assert set(mapping.keys()) == {"7732-18-5", "74-82-8"}
        assert failed == []

    @patch("pubchem.pubchem.safe_request")
    def test_failed_cas(self, mock_req):
        mock_req.return_value = _json_bytes(
            {"Fault": {"Code": "PUGREST.NotFound"}}
        )
        mapping, failed = cas_to_pubchem("NOTEXIST", substance=False)
        assert mapping == {}
        assert "NOTEXIST" in failed

    def test_none_input(self):
        mapping, failed = cas_to_pubchem(None, substance=False)  # type: ignore[arg-type]
        assert mapping == {}
        assert failed == []

    @patch("pubchem.pubchem.safe_request")
    def test_cas_to_cid_wrapper(self, mock_req):
        mock_req.return_value = _json_bytes(
            {"IdentifierList": {"CID": [962]}}
        )
        mapping, failed = cas_to_cid("7732-18-5")
        assert "7732-18-5" in mapping

    @patch("pubchem.pubchem.safe_request")
    def test_cas_to_sid_wrapper(self, mock_req):
        mock_req.return_value = _json_bytes(
            {"IdentifierList": {"SID": [12345]}}
        )
        mapping, failed = cas_to_sid("7732-18-5")
        assert "7732-18-5" in mapping


# ---------------------------------------------------------------------------
# inchikey_to_pubchem
# ---------------------------------------------------------------------------


class TestInchikeyToPubchem:
    @patch("pubchem.pubchem.safe_request")
    def test_returns_cids(self, mock_req):
        mock_req.return_value = _json_bytes(
            {"IdentifierList": {"CID": [962]}}
        )
        result = inchikey_to_pubchem("XLYOFNOQVPJJNP-UHFFFAOYSA-N")
        assert result == [962]

    @patch("pubchem.pubchem.safe_request")
    def test_fault_returns_none(self, mock_req):
        mock_req.return_value = _json_bytes(
            {"Fault": {"Code": "PUGREST.NotFound"}}
        )
        assert inchikey_to_pubchem("AAAAAAAAAAAAAAAA-UHFFFAOYSA-N") is None


# ---------------------------------------------------------------------------
# SMILES_to_pubchem
# ---------------------------------------------------------------------------


class TestSmilesToPubchem:
    @patch("pubchem.pubchem.safe_request")
    def test_simple_smiles(self, mock_req):
        mock_req.return_value = _json_bytes(
            {"IdentifierList": {"CID": [962]}}
        )
        result = SMILES_to_pubchem("O")
        assert result == [962]

    @patch("pubchem.pubchem.safe_request")
    def test_hash_encoding(self, mock_req):
        """The # character should be encoded as %23 in the URL."""
        mock_req.return_value = _json_bytes({"IdentifierList": {"CID": [1]}})
        SMILES_to_pubchem("CC#N")
        called_url = mock_req.call_args[0][0]
        assert "%23" in called_url
        assert "#" not in called_url

    @patch("pubchem.pubchem.safe_request")
    def test_fault_returns_none(self, mock_req):
        mock_req.return_value = _json_bytes({"Fault": {"Code": "PUGREST.NotFound"}})
        assert SMILES_to_pubchem("INVALID") is None


# ---------------------------------------------------------------------------
# get_from_cids
# ---------------------------------------------------------------------------


class TestGetFromCids:
    @patch("pubchem.pubchem.safe_request")
    def test_single_cid(self, mock_req):
        payload = {"PropertyTable": {"Properties": [{"CID": 962, "InChI": "InChI=1S/H2O/h1H2"}]}}
        mock_req.return_value = _json_bytes(payload)
        result = get_from_cids(962, target="property/InChI")
        assert result is not None
        assert result["PropertyTable"]["Properties"][0]["InChI"] == "InChI=1S/H2O/h1H2"

    @patch("pubchem.pubchem.safe_request")
    def test_invalid_json_returns_none(self, mock_req):
        mock_req.return_value = b"not json"
        result = get_from_cids([1, 2], target="synonyms")
        assert result is None


# ---------------------------------------------------------------------------
# get_cids_from_sids
# ---------------------------------------------------------------------------


class TestGetCidsFromSids:
    @patch("pubchem.pubchem.safe_request")
    def test_single_sid(self, mock_req):
        payload = {
            "InformationList": {
                "Information": [{"SID": 12345, "CID": 962}]
            }
        }
        mock_req.return_value = _json_bytes(payload)
        result = get_cids_from_sids(12345)
        assert result == {12345: 962}

    @patch("pubchem.pubchem.safe_request")
    def test_multiple_sids(self, mock_req):
        payload = {
            "InformationList": {
                "Information": [
                    {"SID": 1, "CID": 10},
                    {"SID": 2, "CID": 20},
                ]
            }
        }
        mock_req.return_value = _json_bytes(payload)
        result = get_cids_from_sids([1, 2])
        assert result == {1: 10, 2: 20}


# ---------------------------------------------------------------------------
# cas_to_inchi
# ---------------------------------------------------------------------------


class TestCasToInchi:
    @patch("pubchem.pubchem.safe_request")
    def test_basic_conversion(self, mock_req):
        cid_response = _json_bytes({"IdentifierList": {"CID": [962]}})
        inchi_response = _json_bytes({
            "PropertyTable": {
                "Properties": [{"CID": 962, "InChI": "InChI=1S/H2O/h1H2"}]
            }
        })
        mock_req.side_effect = [cid_response, inchi_response]
        result, failed = cas_to_inchi("7732-18-5")
        assert "7732-18-5" in result
        assert result["7732-18-5"] == "InChI=1S/H2O/h1H2"
        assert failed == []


# ---------------------------------------------------------------------------
# cids_to_cas_and_einecs_and_dtx
# ---------------------------------------------------------------------------


class TestCidsToCasAndEinecs:
    @patch("pubchem.pubchem.safe_request")
    def test_extracts_cas(self, mock_req):
        payload = {
            "InformationList": {
                "Information": [
                    {
                        "CID": 962,
                        "Synonym": ["water", "7732-18-5", "231-791-2"],
                    }
                ]
            }
        }
        mock_req.return_value = _json_bytes(payload)
        result = cids_to_cas_and_einecs_and_dtx([962])
        assert result.get(962, {}).get("CAS") == "7732-18-5"


# ---------------------------------------------------------------------------
# PubchemInputError
# ---------------------------------------------------------------------------


def test_pubchem_input_error_is_attribute_error():
    with pytest.raises(PubchemInputError):
        raise PubchemInputError("test error")

    with pytest.raises(AttributeError):
        raise PubchemInputError("test error")
