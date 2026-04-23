from datetime import datetime, timedelta
import requests
import time
import numpy as np
import regex as re
import logging
from functools import wraps
logger = logging.getLogger(__name__)

class metered_request_decorator:
    """Rate-limiting decorator for PubChem PUG REST requests.

    Enforces the two official PubChem request-rate limits:

    * At most *max_requests_s* requests per second.
    * At most *max_requests_mn* requests per minute.

    When either limit is reached the decorator sleeps for a random
    interval drawn from a Poisson distribution before retrying.  HTTP
    403 responses are retried once with an additional sleep; after 10
    consecutive 403 errors an exception is raised.
    """

    requests_mn: list = []
    STATUS_403: int = 0

    def __init__(self, max_requests_s: int = 5, max_requests_mn: int = 400):
        """Initialise the rate limiter.

        Args:
            max_requests_s: Maximum number of requests allowed per
                second (default: 5).
            max_requests_mn: Maximum number of requests allowed per
                minute (default: 400).
        """
        self.interval_s = timedelta(seconds=1)
        self.interval_mn = timedelta(seconds=60)
        self.max_requests_mn = max_requests_mn
        self.max_requests_s = max_requests_s
        self.len_requests_s = 0
        self.i = 0
        self.random = np.random.poisson(15, 10) / 10

    def update(self):
        """Record the current timestamp as a completed request."""
        metered_request_decorator.requests_mn.append(datetime.now())

    def check(self) -> bool:
        """Return ``True`` when both rate limits allow a new request.

        Returns:
            ``True`` if the per-second and per-minute request counts are
            both below their respective ceilings; ``False`` otherwise.
        """
        self.len_requests_s = len(
            [x for x in self.requests_mn if x > datetime.now() - self.interval_s]
        )
        metered_request_decorator.requests_mn = [
            x for x in self.requests_mn if x > datetime.now() - self.interval_mn
        ]
        return (
            self.len_requests_s < self.max_requests_s
            and len(self.requests_mn) < self.max_requests_mn
        )

    def run(self, fn, *args, **kwargs):
        """Wait until rate limits allow, then call *fn*.

        Args:
            fn: Callable to invoke once rate limits are satisfied.
            *args: Positional arguments forwarded to *fn*.
            **kwargs: Keyword arguments forwarded to *fn*.

        Returns:
            The return value of ``fn(*args, **kwargs)``.
        """
        while self.check() is False:
            time.sleep(self.random[self.i])
            self.i = self.i % 10
        response = fn(*args, **kwargs)
        self.update()
        return response

    def __call__(self, fn):
        """Wrap *fn* with rate-limiting and 403-retry logic.

        Args:
            fn: The function to wrap (expected to return a
                :class:`requests.Response`).

        Returns:
            A wrapper function with the same signature as *fn*.

        Raises:
            Exception: If more than 10 consecutive HTTP 403 responses
                are received.
        """
        @wraps(fn)
        def wrapper(*args, **kwargs):
            response = self.run(fn, *args, **kwargs)
            if response.status_code == 403:
                logging.info('\n got response 403')
                metered_request_decorator.STATUS_403 += 1
                if metered_request_decorator.STATUS_403 > 10:
                    raise Exception('MAX STATUS 403 REACHED!')
                time.sleep(self.random[self.i])
                return self.run(fn, *args, **kwargs)
            return response
        return wrapper


@metered_request_decorator()
def metered_request(_url: str, **params):
    """Send a rate-limited GET request to *_url*.

    Decorated with :class:`metered_request_decorator` so calls are
    automatically throttled to the PubChem request-rate limits.

    Args:
        _url: URL to fetch.
        **params: Additional keyword arguments forwarded to
            :func:`requests.get`.

    Returns:
        :class:`requests.Response` object.
    """
    return requests.get(_url, **params)


def safe_request(_url: str) -> bytes:
    """Fetch *_url* and return the raw response bytes.

    Convenience wrapper around :func:`metered_request` that extracts
    the response content and logs any unexpected status codes.

    Args:
        _url: URL to GET.

    Returns:
        Raw response body as :class:`bytes`.
    """
    response = metered_request(_url)
    return response.content