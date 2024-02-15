from datetime import datetime, timedelta
import requests
import time
import numpy as np
import regex as re
import logging
from functools import wraps
logger = logging.getLogger(__name__)

class metered_request_decorator:
    """
    Prevent sending too many queries to pubchem pug rest

    :params max_q_mn: Maximal number of queries per minute
    :params max_q_s: Maximal number of queries per second
    """
    requests_mn = []
    STATUS_403 = 0
    def __init__(self, max_requests_s:int=5, max_requests_mn:int = 400):
        """
        __init__
        """
        self.interval_s = timedelta(seconds = 1)
        self.interval_mn = timedelta(seconds = 60)
        self.max_requests_mn = max_requests_mn
        self.max_requests_s = max_requests_s
        self.len_requests_s = 0
        self.i = 0
        self.random = np.random.poisson(15,10)/10
    def update(self):
        metered_request_decorator.requests_mn.append(datetime.now())
        #logger.info(f"{datetime.now():%H:%M:%S}: {len(self.requests_mn)} requests")
    def check(self):
        self.len_requests_s = len([x for x in self.requests_mn if x > datetime.now()-self.interval_s])
        metered_request_decorator.requests_mn = [x for x in self.requests_mn if x > datetime.now()-self.interval_mn]
        return self.len_requests_s<self.max_requests_s and len(self.requests_mn)<self.max_requests_mn
    def run(self, fn,*args, **kwargs):
        while self.check() is False:
                time.sleep(self.random[self.i])
                self.i = self.i % 10
        response = fn(*args, **kwargs)
        self.update()
        return response
    def __call__(self, fn):
        def wrapper(*args, **kwargs):
            response = self.run(fn,*args,**kwargs)
            if response.status_code == 403:
                logging.info('\n got response 403')
                metered_request_decorator.STATUS_403+=1
                if metered_request_decorator.STATUS_403 > 10:
                    raise Exception('MAX STATUS 403 REACHED!')
                time.sleep(self.random[self.i])
                return self.run(fn,*args,**kwargs)
            return response
        return wrapper

@metered_request_decorator()
def metered_request(_url:str,**params):
    """
    Call url with limit defined in metered_requests_decorator

    :params _url: url to get
    :params **params: parameters of the request
    """
    return requests.get(_url, **params)

def safe_request(_url:str)->str:
    """
    Return content of the url request, log error messages
    
    :params _url: GET request
    """
    response = metered_request(_url)
    return response._content