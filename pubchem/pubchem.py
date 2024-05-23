from typing import Union
import json
import regex as re
import os
from io import BytesIO
from tqdm import tqdm
import logging
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

def cas_to_pubchem(cas:Union[str,list], substance:bool):
    """
    wrapper for single_cas_to_pubchem

    :params cas: `<str>` or `<list>` of length 3 with cas number to convert or list of cas
    :params substance: `<bool>` True if substance (sid), False if compound (cid)
    """
    fcas = []
    if isinstance(cas, list) and not(len(cas)==3 and len(str(cas[0]))<8\
          and len(str(cas[1]))==2 and len(str(cas[2]))==1):
        #there are multiple cas
        fcas = [format_cas(c) for c in cas]
    elif isinstance(cas,str):
        fcas = [cas]
    elif isinstance(cas,list) and len(cas)==3 and len(str(cas[0]))<8\
          and len(str(cas[1]))==2 and len(str(cas[2]))==1:
        fcas =[format_cas(cas)]
    elif cas is None:
        return {}, []
    # Divide the process into 10 increments wait 1s between each
    return_dict = {}
    failed = []
    for i,_cas in tqdm(enumerate(fcas), total=len(fcas)):
        x=single_cas_to_pubchem(_cas,substance=substance)
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

def cids_to_cas_and_einecs(cids:list, max_query:int = 100):
    """
    Searches for cas and einecs corresponding to cids

    :params cids: list of cids for which to search for cas and einecs
    :params max_query: chunk size of cids sent to Pubchem
    """
    PATcas = r'^(\d){2,7}[-](\d){2}[-](\d)$'
    PATeinecs = r'^(\d){3}[-](\d){3}[-](\d)$'
    cids_chunks = get_chunks(cids, max_query)
    cids_cas_einecs = {}
    for chunk in tqdm(cids_chunks):
        synonyms = get_from_cids(chunk, target = 'synonyms')
        if synonyms is not None:
            for v in synonyms.get('InformationList',{}).get('Information',[]):
                for candidate in v.get('Synonym',[]):
                    cas = re.match(PATcas,candidate)
                    einecs = re.match(PATeinecs,candidate)
                    if cas is not None:
                        cids_cas_einecs.setdefault(v['CID'],{})['CAS'] = cas[0]
                    if einecs is not None:
                        cids_cas_einecs.setdefault(v['CID'],{})['EINECS'] = einecs[0]
    return cids_cas_einecs


def single_cas_to_pubchem(cas,substance:bool=False):
    """Use PubChem REST to convert CAS to CID

    :params cas: `<str>` or `<list>` of length 3 with cas number to convert
    :params substance: `<bool>` True if substance (sid), False if compound (cid)

    :return: cid or sid
    """
    if isinstance(cas, list) and len(cas)==3 and len(str(cas[0]))<8\
          and len(str(cas[1]))==2 and len(str(cas[2]))==1:
        cas = f"{cas[0]}-{cas[1]}-{cas[2]}"
    elif isinstance(cas,str):
        pass
    else:
        logger.info(f"cas_to_cid: {cas} has not a valid format, returning None")
        return None
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
        return data.get('IdentifierList',{}).get('SID',None)
    else:
        return data.get('IdentifierList',{}).get('CID',None)



def get_from_cids(cids, target= None, format = "JSON"):
    """
    Fetch data from PubChem for the cids provided

    :params cids: `<list>` of int or `<int>`
    :params target: `<str>`, e.g. synonyms
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

def cas_to_mols(cas:Union[list,str],cas_cids=None, cas_sids=None, save=None)->dict:
    """
    :params cas_cids: if provided, will skip fetching cids
    :params cas_sids: if provided, will skip fetching sids
    :save: saving path (without extension)

    Returns a dictionary with Chem.MolSupplier for each cas given
    """
    if cas_cids is None:
        cas_cids, failed = cas_to_pubchem(cas, substance = False)
        if save is not None:
            with open(save+"_cas_cids.json","w") as f:
                json.dump(cas_cids,f)
    else:
        failed = cas
    if cas_sids is None:
        cas_sids, failed = cas_to_pubchem(failed, substance = True)
        if save is not None:
            with open(save+"_cas_sids.json","w") as f:
                json.dump(cas_sids,f)
    cids_from_sids = get_cids_from_sids([y for x in cas_sids.values() for y in x])
    reverse_lookup = {y:cas for cas,sids in cas_sids.items() for y in sids}
    for sid,cids in cids_from_sids.items():
        if cids is not None:
            cas = reverse_lookup[sid]
            cas_cids[cas] = cas_cids.setdefault(cas,[])+cids
        else:
            print(f"No cid for {sid}, {reverse_lookup[sid]}")
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
