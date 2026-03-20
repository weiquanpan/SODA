"""Public package interface for SODA."""

from soda._version import __version__
from soda.model import SODA, SODAResult
from soda.scanpy import integrate_adata

__all__ = ["SODA", "SODAResult", "integrate_adata", "__version__"]

