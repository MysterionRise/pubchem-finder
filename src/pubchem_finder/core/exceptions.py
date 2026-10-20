"""Custom exceptions for PubChem Finder."""


class PubChemFinderError(Exception):
    """Base exception for all PubChem Finder errors."""


class MoleculeValidationError(PubChemFinderError):
    """Raised when molecule validation fails."""


class FingerprintGenerationError(PubChemFinderError):
    """Raised when a fingerprint cannot be generated."""


class IndexNotBuilt(PubChemFinderError):
    """No finalized local similarity index exists at the configured location."""


class IndexBuildError(PubChemFinderError):
    """Building or merging the local similarity index failed."""


class DownloadError(PubChemFinderError):
    """Downloading a PubChem bulk file failed or its checksum did not match."""
