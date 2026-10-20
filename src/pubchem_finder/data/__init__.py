"""Data access: PubChem web services (and, later, bulk downloads)."""

from .pug_rest import (
    BadQuery,
    CompoundNotFound,
    Namespace,
    PubChemClient,
    PubChemError,
    PubChemUnavailable,
)

__all__ = [
    "BadQuery",
    "CompoundNotFound",
    "Namespace",
    "PubChemClient",
    "PubChemError",
    "PubChemUnavailable",
]
