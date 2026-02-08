from typing import Union
import json
import regex as re
import os
from io import BytesIO
from tqdm import tqdm
import logging
import pandas as pd
from rdkit import Chem
logger = logging.getLogger(__name__)
from .throttle import safe_request
STATUS_403 = 0
class PubchemInputError(AttributeError):
    pass

def format_cas(_cas:Union[str,list]):
    """Ensure a cas is formatted as 'xxxxxxxx-yy-z"""
    if isinstance(_cas, list) and len(_cas)==3 and len(str(_cas[0]))<8\
        and len(str(_cas[1]))==2 and len(str(_cas[2]))==1:
        return f"{_cas[0]}-{_cas[1]}-{_cas[2]}"
    elif isinstance(_cas,str) and '-' not in _cas:
        xxxx = _cas[:-3]
        z = _cas[-1]
        yy = _cas[-3:-1]
        if len(xxxx)>8:
            logger.error(f'Error with formatting of cas: {_cas}')
        return f'{xxxx}-{yy}-{z}'
    elif isinstance(_cas,str):
        return _cas
    logger.error(f'Error with formatting of cas : {_cas}')

def cas_to_cid(cas:Union[str,list]):
    """
    wrapper for cas_to_pubchem

    :params cas: `<str>` or `<list>` of length 3 with cas number to convert or list of cas
    """
    return cas_to_pubchem(cas,substance=False)

def cas_to_sid(cas:Union[str,list]):
    """
    wrapper for cas_to_pubchem

    :params cas: `<str>` or `<list>` of length 3 with cas number to convert or list of cas
    """
    return cas_to_pubchem(cas,substance=True)  

def inchikey_to_pubchem(inchikey:str):
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/{inchikey}/cids/JSON"
    response = safe_request(url)
    try:
        data = json.loads(response)
    except json.JSONDecodeError:
        logging.info(response)
        print(f'Got error for {inchikey} in json response')
        return None
    if 'Fault' in data.keys():
        logging.info(f'CAS {inchikey} got response {data.get("Fault","")}')
        return None
    return data.get('IdentifierList',{}).get('CID',None)

def cas_to_pubchem(cas:Union[str,list], substance:bool):
    """
    wrapper for single_synonym_to_pubchem

    :params cas: `<str>` or `<list>` of length 3 with cas number to convert or list of cas
    :params substance: `<bool>` True if substance (sid), False if compound (cid)
    """
    fcas = []
    # if isinstance(cas, list) and not(len(cas)==3 and len(str(cas[0]))<8\
    #       and len(str(cas[1]))==2 and len(str(cas[2]))==1):
    #     #there are multiple cas
    #     fcas = [format_cas(c) for c in cas]
    # elif 
    if isinstance(cas,str):
        fcas = [cas]
    elif cas is None:
        return {}, []
    else:
        fcas = cas
    # Divide the process into 10 increments wait 1s between each
    return_dict = {}
    failed = []
    for i,_cas in tqdm(enumerate(fcas), total=len(fcas)):
        x=single_synonym_to_pubchem(_cas,substance=substance)
        if x is None:
            failed.append(_cas)
        else:
            return_dict[_cas] = x
    return return_dict, failed
    
def get_chunks(l:list, n:int):
    """
    Divide list into chunks of maximal size n

    :params l: list
    :params n: max size of each chunk
    """
    try:
        return [l[i:i+n] for i in range(0, len(l), max(1,n))]
    except TypeError as e:
        logging.info(e)
        return []

def cas_to_inchi(cas:Union[str,list], max_query:int = 100):
    """
    Convert CAS number(s) to InChI string(s)

    :params cas: `<str>` or `<list>` of CAS numbers to convert
    :params max_query: chunk size for batch queries
    :return: dictionary mapping CAS to SMILES, and list of failed CAS
    """
    # First try to get CIDs directly from CAS
    cas_cids, failed = cas_to_pubchem(cas, substance=False)
    
    # the next part is not necessary anymore since CIDs are fetched directly from CAS even for substances
    # # For failed CAS, try to get SIDs
    # if failed:
    #     cas_sids, still_failed = cas_to_pubchem(failed, substance=True)
        
    #     # Convert SIDs to CIDs
    #     if cas_sids:
    #         all_sids = [sid for sids in cas_sids.values() for sid in sids]
    #         sids_to_cids_map = get_cids_from_sids(all_sids)
            
    #         # Map back to CAS
    #         for cas_num, sids in cas_sids.items():
    #             cids_from_sids = []
    #             for sid in sids:
    #                 cid = sids_to_cids_map.get(sid)
    #                 if cid:
    #                     if isinstance(cid, list):
    #                         cids_from_sids.extend(cid)
    #                     else:
    #                         cids_from_sids.append(cid)
    #             if cids_from_sids:
    #                 cas_cids[cas_num] = cids_from_sids
    #                 still_failed.remove(cas_num) if cas_num in still_failed else None
    #     failed = still_failed
    
    # Get SMILES for all CIDs
    cas_inchi = {}
    all_cids = list(set([cid for cids in cas_cids.values() for cid in cids]))
    
    if all_cids:
        cids_chunks = get_chunks(all_cids, max_query)
        cid_to_inchi = {}
        
        for chunk in cids_chunks:
            properties = get_from_cids(chunk, target='property/InChI')
            if properties:
                for prop in properties.get('PropertyTable', {}).get('Properties', []):
                    cid = prop.get('CID')
                    inchi = prop.get('InChI')
                    if cid and inchi:
                        cid_to_inchi[cid] = inchi

        # Map InChI back to CAS
        for cas_num, cids in cas_cids.items():
            inchi_list = [cid_to_inchi.get(cid) for cid in cids if cid in cid_to_inchi]
            if inchi_list:
                cas_inchi[cas_num] = inchi_list[0] if len(inchi_list) == 1 else inchi_list

    return cas_inchi, failed

def cids_to_cas_and_einecs(cids:list, max_query:int = 100):
    """Legacy function, replaced by cids_to_cas_and_einecs_and_dtx"""
    return cids_to_cas_and_einecs_and_dtx(cids, max_query)

def cids_to_cas_and_einecs_and_dtx(cids:list, max_query:int = 100):
    """
    Searches for cas and einecs corresponding to cids

    :params cids: list of cids for which to search for cas and einecs
    :params max_query: chunk size of cids sent to Pubchem
    """
    PATcas = r'^(\d){2,7}[-](\d){2}[-](\d)$'
    PATeinecs = r'^(\d){3}[-](\d){3}[-](\d)$'
    PATdtx = r'((?<=DTXSID)\d{5,10})'
    cids_chunks = get_chunks(cids, max_query)
    cids_cas_einecs = {}
    for chunk in tqdm(cids_chunks):
        synonyms = get_from_cids(chunk, target = 'synonyms')
        if synonyms is not None:
            for v in synonyms.get('InformationList',{}).get('Information',[]):
                for candidate in v.get('Synonym',[]):
                    cas = re.findall(PATcas,candidate)
                    einecs = re.findall(PATeinecs,candidate)
                    dtx = re.findall(PATdtx, candidate)
                    if cas is not None and len(cas)>0:
                        cids_cas_einecs.setdefault(v['CID'],{})['CAS'] = cas[0]
                    if einecs is not None and len(einecs)>0:
                        cids_cas_einecs.setdefault(v['CID'],{})['EINECS'] = einecs[0]
                    if dtx is not None and len(dtx)>0:
                        cids_cas_einecs.setdefault(v['CID'],{})['DTXSID'] = dtx[0]
    return cids_cas_einecs

def single_cas_to_pubchem(cas, substance:bool=False):
    """Legacy method"""
    return single_synonym_to_pubchem(cas, substance)

def single_synonym_to_pubchem(cas,substance:bool=False):
    """Use PubChem REST to convert CAS or other synonyms to CID 

    :params cas: `<str>` or `<list>`
    :params substance: `<bool>` True if substance (sid), False if compound (cid)
    :return: cid or sid
    """
    if substance is True:
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/substance/name/{cas}/sids/JSON"
    else:
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{cas}/cids/JSON"
    response = safe_request(url)
    try:
        data = json.loads(response)
    except json.JSONDecodeError:
        logging.info(response)
        print(f'Got error for {cas} in json response')
        return None
    if 'Fault' in data.keys():
        logging.info(f'CAS {cas} got response {data.get("Fault","")}')
        return None
    if substance is True:
        return data.get('IdentifierList',{}).get('SID',[])
    else:
        return data.get('IdentifierList',{}).get('CID',[])


def SMILES_to_pubchem(smiles):
    """Use PubChem REST to convert CAS or other synonyms to CID 

    :params cas: `<str>` or `<list>`
    :params substance: `<bool>` True if substance (sid), False if compound (cid)

    :return: cid or sid
    """
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/{smiles.replace('#','%23')}/cids/JSON"
    response = safe_request(url)
    try:
        data = json.loads(response)
    except json.JSONDecodeError:
        logging.info(response)
        print(f'Got error for {smiles} in json response')
        return None
    if 'Fault' in data.keys():
        logging.info(f'CAS {smiles} got response {data.get("Fault","")}')
        return None
    return data.get('IdentifierList',{}).get('CID',[])


def get_from_cids(cids, target= None, format = "JSON"):
    """
    Fetch data from PubChem for the cids provided

    :params cids: `<list>` of int or `<int>`
    :params target: `<str>`, e.g. synonyms or e.g. property/InChI
    """
    if cids is not None and not isinstance(cids,list):
        cids = [cids]
    fcid = ','.join([str(x) for x in cids])
    if target is not None:
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{fcid}/{target}/{format}"
    else:
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{fcid}/{format}"
    response = safe_request(url)
    try:
        return json.loads(response)
    except Exception as e:
        logger.error(f'Could not fetch {url}, got \n {e}')
        return None
    
def get_mols_from_cas(cas:Union[list,str],substance = False)-> Chem.SDMolSupplier:
    cids, success = cas_to_pubchem(cas,substance = substance)
    return get_mols_from_cids(cids.values())

def get_cids_from_sids(sids:Union[list,int], max_query = 50):
    if sids is not None and not isinstance(sids,list):
        sids = [sids]
    sids_chunks = get_chunks(sids,max_query)
    ret = {}
    for csids in sids_chunks:
        fsid = ','.join([str(x) for x in csids])
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/substance/sid/{fsid}/cids/JSON"
        data = safe_request(url)
        try:
            data = json.loads(data.decode('utf8'))
        except Exception as e:
            print(data)
            raise e
        data = data.get('InformationList',{}).get('Information',[])
        ret.update({x['SID']:x.get('CID',None) for x in data})
    return ret


def get_mols_from_cids(cids:Union[list,int], index = 0, max_query = 500, filename = None) -> Chem.SDMolSupplier:
    """
    Fetch SDF from PubChem for the cids provided

    :params cids: `<list>` of int or `<int>`
    :return: instance of rdkit.Chem  SDMolSupplier
    """
    if cids is not None and not isinstance(cids,list):
        cids = [cids]
    if filename is None:
        filename = f"sdffrompubchem_temp_{index}.sdf"
    cids_chunks = get_chunks(cids,max_query)
    for ccids in cids_chunks:
        fcid = ','.join([str(x) for x in ccids])
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{fcid}/SDF"
        data = safe_request(url)
        if data[0:6]=='<?xml':
            logger.error('in get_mols_from_cids: received an error from server')
            raise PubchemInputError('in get_mols_from_cids: received an error from server')
        with open(filename,'ab') as f:
                f.write(data)
    try:
        suppl = Chem.SDMolSupplier(filename,
                                    sanitize = True,
                                    removeHs = False)
        # suppl = Chem.ForwardSDMolSupplier(stream,
        #                             sanitize = False,
        #                             removeHs = True)
        return suppl, filename
    except Exception as e:
        logger.info(e)
        print(f"Got error {e.__class__} for {index}")
        return [], ''

def get_mols_from_sids(sids:Union[list,int], index = 0, max_query = 200, filename = None) -> Chem.SDMolSupplier:
    """
    Fetch SDF from PubChem for the cids provided

    :params sids: `<list>` of int or `<int>`
    :return: instance of rdkit.Chem  SDMolSupplier
    """
    if sids is not None and not isinstance(sids,list):
        sids = [sids]
    if filename is None:
        filename = f"sdffrompubchem_temp_{index}.sdf"
    sids_chunks = get_chunks(sids,max_query)
    for ssids in sids_chunks:
        fsid = ','.join([str(x) for x in ssids])
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/substance/sid/{fsid}/SDF"
        data = safe_request(url)
        if data[0:6]=='<?xml':
            logger.error('in get_mols_from_sids: received an error from server')
            raise PubchemInputError('in get_mols_from_sids: received an error from server')
        with open(filename,'ab') as f:
                f.write(data)
    try:
        suppl = Chem.SDMolSupplier(filename,
                                    sanitize = True,
                                    removeHs = False)
        # suppl = Chem.ForwardSDMolSupplier(stream,
        #                             sanitize = False,
        #                             removeHs = True)
        return suppl, filename
    except Exception as e:
        logger.info(e)
        print(f"Got error {e.__class__} for {index}")
        return [], ''

def cids_to_mol(cids, filename = None, max_query = 500):
    suppl, _ = get_mols_from_cids(cids, filename = filename, max_query = max_query)
    for i,mol in enumerate(suppl):
        yield mol

def get_cids_from_smiles(smiles):
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/{smiles}/cids/JSON"
    response = safe_request(url)
    data = json.loads(response)
    return data.get('IdentifierList',{}).get('CID',None)

def cas_to_mols(cas:Union[list,str],save=None)->dict:
    """
    :save: saving path (without extension)

    Returns a dictionary with Chem.MolSupplier for each cas given
    """
    cas_cids, failed = cas_to_pubchem(cas, substance = False)
    if save is not None:
        with open(save+"_cas_cids.json","w") as f:
            json.dump(cas_cids,f)
    mols = {}
    files = []
    for i,(cas, cids) in enumerate(cas_cids.items()):
        suppl, filename = get_mols_from_cids(cids, index = i)
        for m in suppl:
            if m is not None:
                m.SetProp('CAS', cas)
                mols.setdefault(cas,[]).append(m) 
        files.append(filename)
    for filename in files:
        try:
            os.remove(filename)
        except FileNotFoundError:
            pass
        except PermissionError:
            pass
    return mols     

def cid_to_cas(cid, get_einecs = False) -> dict:
    """
    Use PubChem REST to convert CID to CAS

    :params cid: `<str>` or `<int>` or `<list>` cid to convert
    :params einec: `<bool>` if True, will return a tupple with cas first and einecs second
    :return: a dictionary of lists with given cids as entries, if only one cid was given
            then the function returns only a list with corresponding cas number as list
    """
    if isinstance(cid,str) or isinstance(cid, int):
        multiple = False
        cid = [cid]
    elif isinstance(cid,list):
        multiple = True
    else:
        logger.info(f"cas_to_cid: {cid} has not a valid format, returning None")
        return None
    data = get_from_cids(cid,'synonyms')
    if data is None:
        return None
    cas = {}
    einecs = {}
    PAT = r"\d{1,10}[-]\d{2}[-]\d"
    EINECS = r'(EINECS\s)?(\d{3}[-]\d{3}[-]\d)'
    for info in data.get("InformationList", {}).get("Information", []):
        icid = info.get('CID',None)
        if icid is not None:
            synonyms = info.get('Synonym',[])
            for synonym in synonyms:
                mat = re.fullmatch(PAT, synonym)
                if mat is not None:
                    cas.setdefault(icid,[]).append(synonym)
                if get_einecs is True:
                    mat = re.fullmatch(EINECS, synonym)
                    if mat is not None:
                        einecs.setdefault(icid,[]).append(mat.group(2))
    if multiple is True and get_einecs is False:
        return cas
    elif get_einecs is False:
        return cas.get(cid[0],[None])
    elif multiple is True:
        return {c: (cas.get(c,[None]),einecs.get(c,[None]))
                    for c in set(cas.keys()).union(set(einecs.keys()))}
    else:
        return {c: (cas.get(c,[None]),einecs.get(c,[None]))
                    for c in set(cas.keys()).union(set(einecs.keys()))}.get(cid[0], None)


def pubchem_pfas_tree(hnid = 5517102, max_tries = 2) -> list:
    """
    Fetch pfas_tree info for specified hnid. To get a list of hid
    for PFAS use the function `getPcHidTree` from `lists.scripts.utility`.


    :params hnid: `<int>` `<str>` with hnid number to process, default to OECD list
    see list lists.data.PubChem_PFAS_Tree_Details.csv
    """
    id_type = 'cids'
    output_format = 'json'
    rest = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/classification/hnid/{str(hnid)}/{id_type}/{output_format}"
    response = safe_request(rest)
    try:
        data = json.loads(response)
    except:
        return None
    if data.get("IdentifierList") is None:
        return None
    return data.get("IdentifierList").get("CID")

def download_cids_by_hnid(hnid:int, folder_path:str, force:bool = False):
    """
    Downloads the list of cids in hnid and saves it in path (last digit of hnid as subfolder)

    :params hnid: hnid of the list
    :params folder_path: the main folder where the json files will be saved
    :params force: whether to force to download all lists (otherwise, lists for which a file already exists will be skiped)
    """
    filename = f"{folder_path}/hnid_{str(hnid).zfill(7)}.json"
    if force is True or not os.path.exists(filename):
        cids_all = pubchem_pfas_tree(hnid=hnid)
        with open(filename, 'w') as f:
            json.dump(cids_all,f)
        return cids_all
    with open(filename, 'r') as f:
            return json.load(f)

def get_HSDB(cids:Union[list,int]):
    """
    returns HSDB id from list of cids    
    """
    data = {}
    PAT = r"(?<=HSDB\s)(\d*)"
    def find_entry(synonyms):
        for synonym in synonyms:
            for match in re.finditer(PAT,synonym):
                return int(match.group())
        return None
    if isinstance(cids,int):
        cids = [cids]
    elif not isinstance(cids,list):
        raise ValueError("cids should be an integer of a list of integers")
    info = get_from_cids(cids,target='synonyms').get("InformationList",{}).get('Information',{})
    for d in info:
        cid = d['CID']
        data[cid] = find_entry(d['Synonym'])
    if len(cids)==1:
        return data[cids[0]]
    return data
    
# Function to process CAS numbers and get SMILES
def cas_to_smiles(cas_list, unique=True, join_smiles = True):
    """
    Convert CAS numbers (or other synonyms, e.g DTXSID) to SMILES, handling multiple molecules per CAS
    and deduplicating by InChI key (if unique =True, else returns list)
    if join_smiles is True, multiple molecules will be converted to SMILES and joined using . (dot).
    """
    all_smiles = []
    processed_cas = {}
    failed_cas = []
    
    for cas in cas_list:
        if pd.isna(cas) or cas == '':
            continue
            
        try:
            print(f"Processing CAS: {cas}")
            mols_dict = cas_to_mols(cas)
            
            if cas in mols_dict:
                mol_list = mols_dict[cas]
                cas_smiles = []
                if unique is True:
                    seen_inchikeys = set()
                    
                    for mol in mol_list:
                        if mol is not None:
                            try:
                                # Generate InChI key for deduplication
                                inchi_key = Chem.MolToInchiKey(mol)
                                
                                if inchi_key not in seen_inchikeys:
                                    smiles = Chem.MolToSmiles(mol)
                                    cas_smiles.append(smiles)
                                    seen_inchikeys.add(inchi_key)
                            except Exception as e:
                                print(f"  Warning: Could not process molecule for CAS {cas}: {e}")
                else:
                    cas_smiles = [Chem.MolToSmiles(m) for m in mol_list]
                if cas_smiles:
                    # Combine multiple SMILES with dots
                    if join_smiles is True:
                        combined_smiles = '.'.join(cas_smiles)
                        processed_cas[cas] = combined_smiles
                        all_smiles.append(combined_smiles)
                    else:
                        processed_cas[cas] = cas_smiles
                    print(f"  Success: {len(cas_smiles)} unique molecules found")
                else:
                    failed_cas.append(cas)
                    processed_cas[cas] = None
                    print(f"  Failed: No valid molecules found")
            else:
                failed_cas.append(cas)
                processed_cas[cas] = None
                print(f"  Failed: CAS not found in PubChem")
                
        except Exception as e:
            print(f"  Error processing CAS {cas}: {e}")
            failed_cas.append(cas)
            processed_cas[cas] = None
    
    return processed_cas, failed_cas
