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
    assert _detailed(PAMOIC_ACID).note is None
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
        assert f"Contains {metal}" in result.note


def test_metal_compounds_drawn_as_ions_are_returned_as_given_too():
    """PubChem, and so this tool's own lookup, writes carboplatin and oxaliplatin
    as separate ions. Only compounds drawn with bonds to the metal were kept:
    these still came back as the ligand, and nitroprusside as hydrogen cyanide."""
    for smiles, metal in (
        ("C1CC(C1)(C(=O)O)C(=O)O.[NH2-].[NH2-].[Pt+2]", "Pt"),  # CID 426756
        ("C1CC[C@H]([C@@H](C1)[NH-])[NH-].C(=O)(C(=O)O)O.[Pt+2]", "Pt"),  # CID 9887053
        ("[C-]#N.[C-]#N.[C-]#N.[C-]#N.[C-]#N.N=O.[Fe+4].[Na+].[Na+]", "Fe"),
        ("[Fe+2].[O-]C(=O)/C=C/C(=O)[O-]", "Fe"),  # ferrous fumarate
        ("[Ag+].Nc1ccc(cc1)S(=O)(=O)[N-]c1ncccn1", "Ag"),  # silver sulfadiazine
    ):
        mol = Chem.MolFromSmiles(smiles)
        result = MoleculeStandardizer().standardize_detailed(mol)
        assert _smiles(result) == Chem.MolToSmiles(mol)
        assert result.removed == []
        assert f"Contains {metal}: without the metal only a ligand" in result.note


def test_group_1_and_2_metals_are_counter_ions_unless_bonded_to_carbon():
    """Sodium acetate drawn with an O-Na bond is a salt; butyllithium is not, and
    was said to have had its lithium "left out" while it was still there."""
    result = _detailed("CC(=O)O[Na]")
    assert (_smiles(result), result.removed, result.note) == (
        "CC(=O)O",
        ["[Na+]"],
        None,
    )
    for smiles, metal in (("[Li]CCCC", "Li"), ("C[Mg]Br", "Mg")):
        mol = Chem.MolFromSmiles(smiles)
        result = MoleculeStandardizer().standardize_detailed(mol)
        assert _smiles(result) == Chem.MolToSmiles(mol)
        assert result.removed == []
        assert result.note.startswith(f"Contains {metal}")


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
    """Lithium carbonate became carbonic acid, with nothing saying the active ion
    was gone."""
    result = _detailed("[Li+].[Li+].[O-]C([O-])=O")
    assert result.removed == ["[Li+]", "[Li+]"]
    assert result.note == (
        "Left out Li: if the metal is the active part, use the input as it is."
    )
    # Sodium, potassium, magnesium and calcium go without comment.
    assert _detailed("[Ca+2].CC(=O)[O-].CC(=O)[O-]").note is None


SALTS = {
    # name: (as given, parent, counter-ion as reported)
    "fosfomycin tromethamine (CID 10199068)": (
        "C[C@H]1[C@H](O1)P(=O)(O)O.C(C(CO)(CO)N)O",
        "C[C@H]1[C@H](O1)P(=O)(O)O",
        "NC(CO)(CO)CO",
    ),
    "diatrizoate meglumine": (
        "CC(=O)Nc1c(I)c(NC(C)=O)c(I)c(C(=O)O)c1I"
        ".CNC[C@H](O)[C@@H](O)[C@H](O)[C@H](O)CO",
        "CC(=O)Nc1c(I)c(NC(C)=O)c(I)c(C(=O)O)c1I",
        "CNC[C@H](O)[C@@H](O)[C@H](O)[C@H](O)CO",
    ),
    "choline salicylate": (
        "C[N+](C)(C)CCO.[O-]C(=O)c1ccccc1O",
        "OC(=O)c1ccccc1O",
        "C[N+](C)(C)CCO",
    ),
    "trolamine salicylate": (
        "OCCN(CCO)CCO.OC(=O)c1ccccc1O",
        "OC(=O)c1ccccc1O",
        "OCCN(CCO)CCO",
    ),
    "lysine acetylsalicylate": (
        "NCCCC[C@H](N)C(=O)O.CC(=O)Oc1ccccc1C(=O)O",
        "CC(=O)Oc1ccccc1C(=O)O",
        "NCCCC[C@H](N)C(=O)O",
    ),
    "cysteamine bitartrate": (
        "NCCS.O[C@H]([C@@H](O)C(=O)O)C(=O)O",
        "NCCS",
        "O[C@H]([C@@H](O)C(=O)O)C(=O)O",
    ),
    "piperazine citrate": (
        "C1CNCCN1.OC(=O)CC(O)(CC(=O)O)C(=O)O",
        "C1CNCCN1",
        "OC(=O)CC(O)(CC(=O)O)C(=O)O",
    ),
    "methenamine hippurate": (
        "C1N2CN3CN1CN(C2)C3.OC(=O)CNC(=O)c1ccccc1",
        "C1N2CN3CN1CN(C2)C3",
        "OC(=O)CNC(=O)c1ccccc1",
    ),
    "dextroamphetamine saccharate": (
        "C[C@H](N)Cc1ccccc1.O[C@@H]([C@H](O)[C@@H](O)C(=O)O)[C@@H](O)C(=O)O",
        "C[C@H](N)Cc1ccccc1",
        "O[C@@H]([C@H](O)[C@@H](O)C(=O)O)[C@@H](O)C(=O)O",
    ),
    "penicillin G benzathine": (
        "CC1(C)S[C@@H]2[C@H](NC(=O)Cc3ccccc3)C(=O)N2[C@H]1C(=O)O"
        ".C(CNCc1ccccc1)NCc1ccccc1",
        "CC1(C)S[C@@H]2[C@H](NC(=O)Cc3ccccc3)C(=O)N2[C@H]1C(=O)O",
        "C(CNCc1ccccc1)NCc1ccccc1",
    ),
}


def test_the_counter_ion_of_a_marketed_salt_is_never_the_parent():
    """RDKit's "largest fragment" counts hydrogens: meglumine (30 atoms) beat
    diatrizoic acid (29) and tromethamine beat fosfomycin, and the counter-ion
    came back as the parent with no note. So did tartrate, citrate, hippurate and
    saccharate next to a small drug."""
    for name, (given, parent, partner) in SALTS.items():
        result = _detailed(given)
        assert _smiles(result) == Chem.CanonSmiles(parent), name
        assert result.removed == [Chem.CanonSmiles(partner)], name
        assert result.note is None, name


def test_another_copy_of_the_parent_was_not_removed():
    """A 2:1 salt: `removed` listed the second salbutamol next to the sulfate, and
    the prompt tells the assistant to report what was removed."""
    salbutamol = "CC(C)(C)NCC(O)c1ccc(O)c(CO)c1"
    result = _detailed(f"{salbutamol}.{salbutamol}.OS(=O)(=O)O")
    assert _smiles(result) == Chem.CanonSmiles(salbutamol)
    assert result.removed == ["O=S(=O)(O)O"]
    acetates = _detailed("[Ca+2].CC(=O)[O-].CC(=O)[O-]")
    assert (_smiles(acetates), acetates.removed) == ("CC(=O)O", ["[Ca+2]"])


def test_an_organic_component_that_is_no_counter_ion_is_pointed_at():
    trimethoprim = "COc1cc(Cc2cnc(N)nc2N)cc(OC)c1OC"
    sulfamethoxazole = "Cc1cc(NS(=O)(=O)c2ccc(N)cc2)no1"
    result = _detailed(f"{trimethoprim}.{sulfamethoxazole}")
    assert result.note == (
        f"Left out {Chem.CanonSmiles(sulfamethoxazole)}, which is not a common "
        "counter-ion or solvent: the input may be a mixture, and the parent returned "
        "is only its largest component."
    )


def test_solvents_and_counter_ions_alone():
    """Pamoic acid with ethanol or DMSO came back as the solvent."""
    for solvent in ("CCO", "CS(C)=O"):
        result = _detailed(f"{solvent}.{PAMOIC_ACID}")
        assert _smiles(result) == Chem.CanonSmiles(PAMOIC_ACID)
        assert result.removed == [solvent]
        assert result.note == (
            "Every organic component is a common counter-ion or solvent; the "
            "largest was kept."
        )
    # One of them on its own, or with its metal, is simply the compound.
    assert _detailed("CCO").note is None
    stearate = "CCCCCCCCCCCCCCCCCC(=O)[O-]"
    assert _detailed(f"{stearate}.{stearate}.[Mg+2]").note is None


def test_a_removed_sulfoxide_is_reported_neutral_too():
    """A DMSO solvate: `removed` must show what was given, not C[S+](C)[O-]."""
    result = _detailed("CS(C)=O.CC(=O)Oc1ccccc1C(=O)O")
    assert _smiles(result) == "CC(=O)Oc1ccccc1C(=O)O"
    assert result.removed == ["CS(C)=O"]
