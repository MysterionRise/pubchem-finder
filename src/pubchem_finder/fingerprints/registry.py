"""Fingerprint registry for dynamic fingerprint selection."""

from ..core.fingerprint import FingerprintGenerator
from .maccs import MACCSFingerprint
from .morgan import MorganFingerprint


class FingerprintRegistry:
    """
    Registry pattern for fingerprint generators.

    This makes it easy to add new fingerprint types without modifying
    existing code. Future ML-based fingerprints can be registered here.

    Example:
        >>> from pubchem_finder.fingerprints import fingerprint_registry
        >>> gen = fingerprint_registry.get("morgan_2")
        >>> gen.name
        'morgan_ecfp_r2_2048'
    """

    def __init__(self) -> None:
        self._generators: dict[str, FingerprintGenerator] = {}
        self._register_defaults()

    def _register_defaults(self) -> None:
        """Register standard fingerprints."""
        # Morgan variants
        self.register("morgan_2", MorganFingerprint(radius=2, n_bits=2048))
        self.register("morgan_3", MorganFingerprint(radius=3, n_bits=2048))
        self.register("ecfp4", MorganFingerprint(radius=2, n_bits=2048))  # Alias
        self.register("ecfp6", MorganFingerprint(radius=3, n_bits=2048))  # Alias

        # Feature-based Morgan (FCFP)
        self.register(
            "fcfp4", MorganFingerprint(radius=2, n_bits=2048, use_features=True)
        )

        # MACCS keys
        self.register("maccs", MACCSFingerprint())

    def register(self, name: str, generator: FingerprintGenerator) -> None:
        """
        Register a new fingerprint generator.

        Args:
            name: Short name for the fingerprint
            generator: FingerprintGenerator instance
        """
        self._generators[name] = generator

    def get(self, name: str) -> FingerprintGenerator:
        """
        Get fingerprint generator by name.

        Args:
            name: Fingerprint name

        Returns:
            FingerprintGenerator instance

        Raises:
            KeyError: If fingerprint name not found
        """
        if name not in self._generators:
            available = ", ".join(self.list_available())
            raise KeyError(f"Unknown fingerprint type: {name}. Available: {available}")
        return self._generators[name]

    def list_available(self) -> list[str]:
        """List all registered fingerprint types."""
        return sorted(self._generators.keys())

    def __contains__(self, name: str) -> bool:
        """Check if fingerprint is registered."""
        return name in self._generators


# Global registry instance
fingerprint_registry = FingerprintRegistry()
