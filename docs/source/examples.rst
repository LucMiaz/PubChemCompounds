Usage Examples
==============

This page shows practical examples for the most common tasks.
All examples assume the package has been installed and imported as::

    import pubchem_compounds as pc


Identifier lookups
------------------

**CAS to CID (single)**

.. code-block:: python

    mapping, failed = pc.cas_to_cid("7732-18-5")  # water
    print(mapping)   # {'7732-18-5': [962]}

**CAS to CID (batch)**

.. code-block:: python

    cas_list = ["7732-18-5", "74-82-8", "71-43-2"]
    mapping, failed = pc.cas_to_cid(cas_list)
    # mapping == {'7732-18-5': [962], '74-82-8': [297], '71-43-2': [241]}
    if failed:
        print("Could not resolve:", failed)

**CAS to SID**

.. code-block:: python

    mapping, failed = pc.cas_to_sid("7732-18-5")
    print(mapping)   # {'7732-18-5': [<SID>, ...]}

**InChIKey to CID**

.. code-block:: python

    cids = pc.inchikey_to_pubchem("XLYOFNOQVPJJNP-UHFFFAOYSA-N")
    print(cids)  # [962]

**SMILES to CID**

.. code-block:: python

    cids = pc.SMILES_to_pubchem("O")   # water
    print(cids)  # [962]

    # Triple bonds are URL-encoded automatically:
    cids = pc.SMILES_to_pubchem("CC#N")  # acetonitrile
    print(cids)

**CID to CAS (reverse lookup)**

.. code-block:: python

    cas = pc.cid_to_cas(962)
    print(cas)  # ['7732-18-5']

    # Extract CAS + EINECS for multiple CIDs:
    result = pc.cids_to_cas_and_einecs_and_dtx([962, 297])
    # {962: {'CAS': '7732-18-5', 'EINECS': '231-791-2'}, ...}


Batch property fetching
-----------------------

.. code-block:: python

    data = pc.get_from_cids(
        [962, 297, 241],
        target="property/MolecularFormula,MolecularWeight"
    )
    for prop in data["PropertyTable"]["Properties"]:
        print(prop["CID"], prop["MolecularFormula"], prop["MolecularWeight"])
    # 962  H2O   18.0
    # 297  CH4   16.0
    # 241  C6H6  78.1

**Fetch synonyms**

.. code-block:: python

    data = pc.get_from_cids([962], target="synonyms")
    synonyms = data["InformationList"]["Information"][0]["Synonym"]
    print(synonyms[:5])


InChI / SMILES conversion
--------------------------

**CAS to InChI**

.. code-block:: python

    result, failed = pc.cas_to_inchi(["7732-18-5", "74-82-8"])
    print(result["7732-18-5"])  # InChI=1S/H2O/h1H2

**CAS to SMILES**

.. code-block:: python

    processed, failed = pc.cas_to_smiles(["7732-18-5", "74-82-8"])
    print(processed["7732-18-5"])  # 'O'
    print(processed["74-82-8"])    # 'C'

**DTXSID to SMILES**

.. code-block:: python

    processed, failed = pc.dtxsid_to_smiles(["DTXSID9020584"])
    print(processed["DTXSID9020584"])

**Any synonym to SMILES** (CAS, DTXSID, compound name, …)

.. code-block:: python

    processed, failed = pc.synonyms_to_smiles(
        ["water", "methane", "benzene"],
        unique=True,        # deduplicate by InChIKey
        join_smiles=True,   # join multi-structure results with '.'
    )


Mol / SDF retrieval (requires RDKit)
-------------------------------------

**CAS to RDKit molecules**

.. code-block:: python

    mols = pc.cas_to_mols(["7732-18-5", "74-82-8"])
    for cas, mol_list in mols.items():
        for mol in mol_list:
            if mol is not None:
                print(cas, mol.GetNumAtoms(), mol.GetNumBonds())

**CIDs to mol list**

.. code-block:: python

    mols, sdf_path = pc.get_mols_from_cids([962, 297, 241])
    for mol in mols:
        if mol is not None:
            print(mol.GetProp("_Name"))

**Yield molecules one at a time (generator)**

.. code-block:: python

    for mol in pc.cids_to_mol([962, 297, 241]):
        if mol is not None:
            print(mol.GetNumAtoms())

**SIDs to molecules**

.. code-block:: python

    mols, _ = pc.get_mols_from_sids([12345, 67890])


PFAS / Classification tree
---------------------------

**Fetch all CIDs from the OECD PFAS list node**

.. code-block:: python

    cids = pc.pubchem_pfas_tree()          # default hnid = 5517102
    print(f"{len(cids)} PFAS CIDs")

**Custom hierarchy node**

.. code-block:: python

    cids = pc.pubchem_pfas_tree(hnid=7453174)

**Cache to disk and re-use**

.. code-block:: python

    import os
    cache_dir = "my_cache"
    os.makedirs(cache_dir, exist_ok=True)

    cids = pc.download_cids_by_hnid(5517102, folder_path=cache_dir)
    # Subsequent calls read from disk unless force=True:
    cids = pc.download_cids_by_hnid(5517102, folder_path=cache_dir)


Rate limiting
-------------

The built-in :func:`~pubchem_compounds.safe_request` honours PubChem's request
limits automatically.  You do not need to add ``time.sleep()`` calls in loops.

.. code-block:: python

    # This loop is safe — the rate limiter throttles requests automatically.
    for cid in my_cid_list:
        data = pc.get_from_cids(cid, target="property/IsomericSMILES")
        smiles = data["PropertyTable"]["Properties"][0]["IsomericSMILES"]
