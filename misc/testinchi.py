#fnord

from ccdc.io import MoleculeReader
rdr = MoleculeReader('CSD')
from inchi import InChIGenerator
gen = InChIGenerator()
mol = rdr.molecule('ZIYZEP')
inchi = gen.generate(rdr[0], include_stereo=True)
print inchi.inchi
