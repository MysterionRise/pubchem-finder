"""Tests for molecule standardizer."""

from rdkit import Chem

from pubchem_finder.chemistry import MoleculeStandardizer


def test_standardizer_removes_salts(ethanol_with_salt_smiles):
    """Test that standardization removes salt counterions."""
    standardizer = MoleculeStandardizer()
    result = standardizer.standardize_smiles(ethanol_with_salt_smiles)
    assert result == "CCO"


def test_standardizer_removes_hydrogens(ethanol_mol):
    """Test that standardization removes explicit hydrogens."""
    # Add explicit hydrogens
    mol_with_h = Chem.AddHs(ethanol_mol)
    assert mol_with_h.GetNumAtoms() > ethanol_mol.GetNumAtoms()

    standardizer = MoleculeStandardizer()
    std_mol = standardizer.standardize(mol_with_h)

    # Should have removed explicit H
    assert std_mol.GetNumAtoms() == ethanol_mol.GetNumAtoms()


def test_standardizer_largest_fragment():
    """Test that standardizer keeps largest fragment."""
    # Ethanol with sodium chloride salt
    smiles = "CCO.[Na+].[Cl-]"
    standardizer = MoleculeStandardizer()
    result = standardizer.standardize_smiles(smiles)
    assert result == "CCO"


def test_standardizer_neutralize_charges():
    """Test charge neutralization."""
    # Carboxylic acid anion
    smiles = "CC([O-])=O"
    standardizer = MoleculeStandardizer(neutralize=True)
    result = standardizer.standardize_smiles(smiles)

    # Should be neutralized to carboxylic acid
    mol = Chem.MolFromSmiles(result)
    total_charge = sum([atom.GetFormalCharge() for atom in mol.GetAtoms()])
    assert total_charge == 0


def test_standardizer_with_neutralize_false():
    """Test standardizer without neutralization."""
    smiles = "CC([O-])=O"
    standardizer = MoleculeStandardizer(neutralize=False)
    result = standardizer.standardize_smiles(smiles)

    mol = Chem.MolFromSmiles(result)
    total_charge = sum([atom.GetFormalCharge() for atom in mol.GetAtoms()])
    # Charge should be preserved
    assert total_charge == -1


def test_standardizer_with_invalid_smiles(invalid_smiles):
    """Test standardizer with invalid SMILES."""
    standardizer = MoleculeStandardizer()
    result = standardizer.standardize_smiles(invalid_smiles)
    assert result is None


def test_standardizer_with_none_mol():
    """Test standardizer with None molecule."""
    standardizer = MoleculeStandardizer()
    result = standardizer.standardize(None)
    assert result is None


def test_standardizer_canonicalization(benzene_smiles):
    """Test that different representations give same canonical form."""
    # Different valid SMILES representations of benzene
    benzene_variants = [
        "c1ccccc1",  # Aromatic notation
        "C1=CC=CC=C1",  # Kekulé form
        "c1cccc2c1ccc2",  # Naphthalene has multiple forms too
    ]

    standardizer = MoleculeStandardizer()

    # Test that aromatic and Kekulé benzene give same result
    benzene_aromatic = standardizer.standardize_smiles(benzene_variants[0])
    benzene_kekule = standardizer.standardize_smiles(benzene_variants[1])

    assert benzene_aromatic is not None
    assert benzene_kekule is not None
    assert benzene_aromatic == benzene_kekule


def test_standardizer_neutralizes_salt_with_charged_parent():
    """The parent is chosen before neutralization, so the counter-ion cannot
    keep the carboxylate charged."""
    standardizer = MoleculeStandardizer()
    assert standardizer.standardize_smiles("CC(=O)[O-].[Na+]") == "CC(=O)O"


def test_standardizer_keeps_charge_of_salt_parent_when_not_neutralizing():
    standardizer = MoleculeStandardizer(neutralize=False)
    assert standardizer.standardize_smiles("CC(=O)[O-].[Na+]") == "CC(=O)[O-]"


def test_standardizer_prefers_organic_fragment():
    """A small organic fragment beats a larger inorganic one."""
    standardizer = MoleculeStandardizer()
    result = standardizer.standardize_smiles("CCO.O=S(=O)(O)O")
    assert result is not None
    assert "S" not in result or result == "CCO"


# --- Found testing the released 2.2.1 against RDKit and PubChem -------------------

PAMOIC_ACID = "O=C(O)c1cc2ccccc2c(Cc2c(O)c(C(=O)O)cc3ccccc23)c1O"
PYRANTEL = "CN1CCCN=C1/C=C/c1cccs1"
OMEPRAZOLE = "COc1ccc2nc(S(=O)Cc3ncc(C)c(OC)c3C)[nH]c2c1"


def _detailed(smiles):
    return MoleculeStandardizer().standardize_detailed(Chem.MolFromSmiles(smiles))


def _smiles(result):
    return Chem.MolToSmiles(result.mol)


def test_sulfoxides_stay_neutral():
    """RDKit's cleanup draws a sulfoxide [S+]-[O-]. The "standardized" omeprazole
    then scored 0.63 against omeprazole in the index (which holds PubChem's S=O),
    found no substructure match, and gained a Brenk alert for a charged sulfur."""
    for smiles in (OMEPRAZOLE, "CS(C)=O", "NC(=O)CS(=O)C(c1ccccc1)c1ccccc1"):
        result = _detailed(smiles)
        assert _smiles(result) == Chem.CanonSmiles(smiles)
        assert all(atom.GetFormalCharge() == 0 for atom in result.mol.GetAtoms())
    # Given charged, it comes back in the same form as everything else.
    charged = "COc1ccc2nc([S+]([O-])Cc3ncc(C)c(OC)c3C)[nH]c2c1"
    assert _smiles(_detailed(charged)) == Chem.CanonSmiles(OMEPRAZOLE)


def test_sulfoxide_stereo_survives():
    esomeprazole = "CC1=CN=C(C(=C1OC)C)C[S@](=O)C2=NC3=C(N2)C=C(C=C3)OC"
    other = esomeprazole.replace("[S@]", "[S@@]")
    assert _smiles(_detailed(esomeprazole)) == Chem.CanonSmiles(esomeprazole)
    assert _smiles(_detailed(other)) == Chem.CanonSmiles(other)
    assert Chem.CanonSmiles(esomeprazole) != Chem.CanonSmiles(other)


def test_sulfones_and_sulfonium_are_left_alone():
    for smiles in ("CS(C)(=O)=O", "CS(=O)(=O)O", "C[S+](C)C", "C[S+](C)(C)=O"):
        assert _smiles(_detailed(smiles)) == Chem.CanonSmiles(smiles)


def test_no_demo_drug_gains_a_charge():
    """15 of the 1,932 demo drugs (the sulfoxides) used to come back with two."""
    import csv
    from importlib.resources import files

    path = files("pubchem_finder.index") / "demo_data" / "fda_approved_drugs.tsv"
    with path.open(encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert len(rows) > 1900
    standardizer = MoleculeStandardizer()

    def charged(mol):
        return sum(1 for atom in mol.GetAtoms() if atom.GetFormalCharge())

    changed = []
    for row in rows:
        mol = Chem.MolFromSmiles(row["smiles"])
        result = standardizer.standardize_detailed(mol)
        assert result is not None, row["cid"]
        assert charged(result.mol) <= charged(mol), row["cid"]
        if Chem.MolToSmiles(result.mol) != Chem.MolToSmiles(mol):
            changed.append(row["cid"])
    # The set holds parent structures, so almost nothing should change.
    assert len(changed) <= 5, changed


def test_pamoate_is_the_counter_ion_not_the_parent():
    """Pamoic acid has 28 heavy atoms, pyrantel 14: "the largest fragment" was the
    counter-ion, and the drug was dropped without a word."""
    result = _detailed(f"{PYRANTEL}.{PAMOIC_ACID}")
    assert _smiles(result) == Chem.CanonSmiles(PYRANTEL)
    assert result.removed == [Chem.CanonSmiles(PAMOIC_ACID)]
    # As a salt of charged components too.
    salt = "CN1CCC[NH+]=C1/C=C/c1cccs1.[O-]C(=O)c1cc2ccccc2c(Cc2c(O)c(C(=O)O)cc3ccccc23)c1O"
    assert _smiles(_detailed(salt)) == Chem.CanonSmiles(PYRANTEL)
    # On its own, or with only its sodium, it is the compound.
    assert _smiles(_detailed(PAMOIC_ACID)) == Chem.CanonSmiles(PAMOIC_ACID)
    disodium = "[Na+].[Na+].[O-]C(=O)c1cc2ccccc2c(Cc2c(O)c(C(=O)[O-])cc3ccccc23)c1O"
    result = _detailed(disodium)
    assert _smiles(result) == Chem.CanonSmiles(PAMOIC_ACID)
    assert result.removed == ["[Na+]", "[Na+]"]


def test_metal_compounds_are_returned_as_given():
    """Cleanup disconnects metals, so cisplatin became "N" (ammonia), carboplatin
    its dicarboxylic acid and oxaliplatin its diamine."""
    for smiles, metal in (
        ("N.N.Cl[Pt]Cl", "Pt"),
        ("C1CC2(C1)C(=O)O[Pt]OC2=O.N.N", "Pt"),
        ("C1CC[C@H]2[C@@H](C1)N[Pt]3(N2)OC(=O)C(=O)O3", "Pt"),
        ("Cl[Zn]Cl", "Zn"),
        ("O=[Ti]=O", "Ti"),
    ):
        mol = Chem.MolFromSmiles(smiles)
        result = MoleculeStandardizer().standardize_detailed(mol)
        assert _smiles(result) == Chem.MolToSmiles(mol)
        assert result.removed == []
        assert f"Contains bonded {metal}" in result.note


def test_inputs_without_an_organic_part_are_returned_as_given():
    """Sodium chloride became HCl, potassium iodide HI, sodium fluoride [Na+]."""
    for smiles in ("[Na+].[Cl-]", "[K+].[I-]", "[Na+].[F-]", "[Ca]", "O"):
        result = _detailed(smiles)
        assert _smiles(result) == Chem.CanonSmiles(smiles)
        assert "No organic component" in result.note


def test_removed_fragments_are_reported():
    assert _detailed("CC(=O)[O-].[Na+]").removed == ["[Na+]"]
    assert _detailed("CC(=O)[O-].[Na+]").note is None
    assert _detailed("CCO").removed == []
    # A two-drug mixture: one component is kept, and the other is named.
    trimethoprim = "COc1cc(Cc2cnc(N)nc2N)cc(OC)c1OC"
    sulfamethoxazole = "Cc1cc(NS(=O)(=O)c2ccc(N)cc2)no1"
    result = _detailed(f"{trimethoprim}.{sulfamethoxazole}")
    assert _smiles(result) == Chem.CanonSmiles(trimethoprim)
    assert result.removed == [Chem.CanonSmiles(sulfamethoxazole)]


def test_a_metal_left_out_is_pointed_at():
    """Ferrous fumarate became fumaric acid and lithium carbonate carbonic acid,
    with nothing saying the active ion was gone."""
    result = _detailed("[Fe+2].[O-]C(=O)/C=C/C(=O)[O-]")
    assert result.removed == ["[Fe+2]"]
    assert "Left out Fe" in result.note
    assert "Left out Li" in _detailed("[Li+].[Li+].[O-]C([O-])=O").note
