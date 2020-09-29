from indigo import *

indigo = Indigo()

m1 = indigo.loadMolecule(
    "O=C1CC(CN1C1=CC=CC=C1F)C(=O)NCC1=CC(=CC=C1)NC(=O)C1CC1"
)
print(m1.canonicalSmiles())
print(m1.fingerprint(type='sim').oneBitsList())
