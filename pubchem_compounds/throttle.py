from datetime import datetime, timedelta
import requests
import time
import numpy as np
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
    403 (request-rate block) and 503 (server busy) responses are retried
    after an additional sleep; after 10 consecutive such responses an
    exception is raised.
    """

    #: Status codes PubChem uses for "slow down" / "busy": retried.
    RETRY_STATUSES = (403, 503)
    MAX_RETRIES = 10

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
            # Cycle through the 10 pre-drawn intervals (the index used to stay at 0).
            self.i = (self.i + 1) % 10
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
            Exception: If more than :attr:`MAX_RETRIES` consecutive HTTP
                403/503 responses are received.
        """
        @wraps(fn)
        def wrapper(*args, **kwargs):
            response = self.run(fn, *args, **kwargs)
            while response.status_code in self.RETRY_STATUSES:
                logger.warning("PubChem returned HTTP %s, retrying", response.status_code)
                metered_request_decorator.STATUS_403 += 1
                if metered_request_decorator.STATUS_403 > self.MAX_RETRIES:
                    raise Exception(
                        f"PubChem returned HTTP {response.status_code} "
                        f"{self.MAX_RETRIES} times in a row")
                time.sleep(self.random[self.i])
                self.i = (self.i + 1) % 10
                response = self.run(fn, *args, **kwargs)
            # The counter is for *consecutive* refusals: a success resets it,
            # otherwise ten refusals spread over a long session would raise.
            metered_request_decorator.STATUS_403 = 0
            return response
        return wrapper


#: Seconds before a request to PubChem/CompTox is abandoned.
REQUEST_TIMEOUT = 60


@metered_request_decorator()
def metered_request(_url: str, **params):
    """Send a rate-limited GET request to *_url*.

    Decorated with :class:`metered_request_decorator` so calls are
    automatically throttled to the PubChem request-rate limits.

    Args:
        _url: URL to fetch.
        **params: Additional keyword arguments forwarded to
            :func:`requests.get`.  ``timeout`` defaults to
            :data:`REQUEST_TIMEOUT` seconds so a stalled connection cannot
            block forever.

    Returns:
        :class:`requests.Response` object.
    """
    params.setdefault("timeout", REQUEST_TIMEOUT)
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