"""Custom exceptions for PubChem Finder."""


class PubChemFinderError(Exception):
    """Base exception for all PubChem Finder errors."""


class MoleculeValidationError(PubChemFinderError):
    """Raised when molecule validation fails."""


class FingerprintGenerationError(PubChemFinderError):
    """Raised when a fingerprint cannot be generated."""
