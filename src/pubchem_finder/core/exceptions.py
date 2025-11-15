"""Custom exceptions for PubChem Finder."""


class PubChemFinderError(Exception):
    """Base exception for all PubChem Finder errors."""

    pass


class MoleculeValidationError(PubChemFinderError):
    """Raised when molecule validation fails."""

    pass


class FingerprintGenerationError(PubChemFinderError):
    """Raised when fingerprint generation fails."""

    pass


class DatabaseError(PubChemFinderError):
    """Raised when database operations fail."""

    pass


class IngestionError(PubChemFinderError):
    """Raised when data ingestion fails."""

    pass
