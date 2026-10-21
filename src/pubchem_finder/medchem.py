"""Medicinal-chemistry triage: drug-likeness, synthesizability, structural alerts.

Everything here runs offline on RDKit:

- QED (Bickerton et al. 2012), a 0-1 drug-likeness score.
- Synthetic accessibility (Ertl & Schuffenhauer 2009) from RDKit's Contrib
  ``SA_Score`` module, 1 (easy) to 10 (hard).
- Lipinski's Rule of Five and Veber's oral-bioavailability rules.
- Structural alerts from the PAINS, Brenk and NIH filter sets, with the
  matched atoms so they can be highlighted.
- The Bemis-Murcko scaffold, also used to group similarity-search hits.
"""

import importlib.util
import logging
import os
import threading
from collections.abc import Sequence
from functools import cache
from types import ModuleType

from rdkit import Chem, RDConfig
from rdkit.Chem import QED, FilterCatalog
from rdkit.Chem.Scaffolds import MurckoScaffold

from .chemistry import calculate_properties
from .models import MedchemProfile, ScaffoldGroup, SimilarityHit, StructuralAlert

logger = logging.getLogger(__name__)

_CATALOGS = (
    FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS,
    FilterCatalog.FilterCatalogParams.FilterCatalogs.BRENK,
    FilterCatalog.FilterCatalogParams.FilterCatalogs.NIH,
)
_lock = threading.Lock()


@cache
def _filter_catalog() -> FilterCatalog.FilterCatalog:
    params = FilterCatalog.FilterCatalogParams()
    for catalog in _CATALOGS:
        params.AddCatalog(catalog)
    return FilterCatalog.FilterCatalog(params)


@cache
def _sascorer() -> ModuleType | None:
    """RDKit's Contrib SA_Score module, or None if this RDKit build lacks it."""
    path = os.path.join(RDConfig.RDContribDir, "SA_Score", "sascorer.py")
    spec = importlib.util.spec_from_file_location("pubchem_finder_sascorer", path)
    if spec is None or spec.loader is None or not os.path.exists(path):
        logger.warning("RDKit Contrib SA_Score not found; sa_score will be None")
        return None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def sa_score(mol: Chem.Mol) -> float | None:
    """Synthetic accessibility score, or None if unavailable."""
    with _lock:  # the module lazily loads its fragment table on first use
        module = _sascorer()
        if module is None:
            return None
        return round(float(module.calculateScore(mol)), 2)


def structural_alerts(mol: Chem.Mol) -> list[StructuralAlert]:
    """PAINS, Brenk and NIH alerts matched by ``mol``."""
    alerts = []
    for entry in _filter_catalog().GetMatches(mol):
        props = set(entry.GetPropList())
        atoms = sorted(
            {
                pair[1]
                for match in entry.GetFilterMatches(mol)
                for pair in match.atomPairs
            }
        )
        alerts.append(
            StructuralAlert(
                catalog=entry.GetProp("FilterSet")
                if "FilterSet" in props
                else "unknown",
                name=entry.GetDescription(),
                atoms=atoms,
                reference=entry.GetProp("Reference") if "Reference" in props else None,
            )
        )
    return alerts


def murcko_scaffold(mol: Chem.Mol) -> str:
    """Bemis-Murcko scaffold SMILES ('' for acyclic molecules)."""
    return str(MurckoScaffold.MurckoScaffoldSmiles(mol=mol))  # type: ignore[no-untyped-call]


def profile(mol: Chem.Mol, input_smiles: str) -> MedchemProfile:
    """Full triage profile for one molecule."""
    props = calculate_properties(mol)
    violations = []
    if props["molecular_weight"] > 500:
        violations.append("molecular weight > 500")
    if props["logp"] > 5:
        violations.append("LogP > 5")
    if props["num_h_donors"] > 5:
        violations.append("H-bond donors > 5")
    if props["num_h_acceptors"] > 10:
        violations.append("H-bond acceptors > 10")
    alerts = structural_alerts(mol)
    return MedchemProfile(
        input_smiles=input_smiles,
        canonical_smiles=Chem.MolToSmiles(mol),
        qed=round(float(QED.qed(mol)), 3),  # type: ignore[no-untyped-call]
        sa_score=sa_score(mol),
        molecular_weight=props["molecular_weight"],
        logp=props["logp"],
        num_h_donors=int(props["num_h_donors"]),
        num_h_acceptors=int(props["num_h_acceptors"]),
        tpsa=props["tpsa"],
        num_rotatable_bonds=int(props["num_rotatable_bonds"]),
        lipinski_violations=violations,
        veber_ok=props["num_rotatable_bonds"] <= 10 and props["tpsa"] <= 140,
        alerts=alerts,
        pains_free=not any(a.catalog.startswith("PAINS") for a in alerts),
        murcko_scaffold=murcko_scaffold(mol),
    )


def group_by_scaffold(hits: Sequence[SimilarityHit]) -> list[ScaffoldGroup]:
    """Group hits by Murcko scaffold, best-scoring group first; order within a
    group follows the input order."""
    groups: dict[str, ScaffoldGroup] = {}
    for hit in hits:
        mol = Chem.MolFromSmiles(hit.smiles) if hit.smiles else None
        key = murcko_scaffold(mol) if mol is not None else "?"
        group = groups.get(key)
        if group is None:
            groups[key] = ScaffoldGroup(
                scaffold=key, cids=[hit.cid], best_similarity=hit.similarity
            )
        else:
            group.cids.append(hit.cid)
            group.best_similarity = max(group.best_similarity, hit.similarity)
    return sorted(groups.values(), key=lambda g: (-g.best_similarity, -len(g.cids)))
