from indigo import *

indigo = Indigo()

m1 = indigo.loadMolecule("OC(Cl)=[C@]=C(C)F")
print(m1.canonicalSmiles())
print(m1.fingerprint(type='sim').oneBitsList())
