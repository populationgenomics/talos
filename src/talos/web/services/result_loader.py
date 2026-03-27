"""
Caching loader for ResultData JSON files.
"""

from collections import OrderedDict

from loguru import logger

from talos.models import ResultData
from talos.utils import read_json_from_path


class ResultCache:
    """LRU cache for parsed ResultData objects, keyed by file path."""

    def __init__(self, max_size: int = 10):
        self._cache: OrderedDict[str, ResultData] = OrderedDict()
        self._max_size = max_size

    def get(self, json_path: str) -> ResultData:
        if json_path in self._cache:
            self._cache.move_to_end(json_path)
            return self._cache[json_path]

        logger.info(f'Loading result file: {json_path}')
        result = read_json_from_path(json_path, return_model=ResultData)
        if result is None:
            raise FileNotFoundError(f'Result file not found or invalid: {json_path}')

        if len(self._cache) >= self._max_size:
            evicted_key, _ = self._cache.popitem(last=False)
            logger.debug(f'Evicted cached result: {evicted_key}')

        self._cache[json_path] = result
        return result

    def invalidate(self, json_path: str):
        self._cache.pop(json_path, None)

    def clear(self):
        self._cache.clear()
