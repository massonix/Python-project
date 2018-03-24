##hemoglobin OK

from Bio.PDB import *
import sys
import os
from Bio import pairwise2

if len(sys.argv) < 2:
    inter = None
else:
    inter = sys.argv[1]

if os.path.isdir(inter):
    files_list = []
    for files in os.listdir(inter):
        if files.endswith(".pdb"):
            files_list.append(os.path.join(inter, files))


def GetStructures(pdbfile):
    """
    Given a pdbfile, gets it's structure.
    :param pdbfile: pdb file to open
    :return: structure object
    """
    parser = PDBParser()
    structure = parser.get_structure(pdbfile[0:-4], pdbfile)
    return structure

class Superimpose(object):

    def __init__(self, structure1, structure2):
        self.structure1 = structure1
        self.structure2 = structure2
        self.si = Superimposer()

    def SuperimposeStructures(self):
        """
        :param structure1: first structure to superimpose
        :param structure2: second structure to superimpose
        :return: rmsd object
            """
        atoms_a = list(self.structure1.get_atoms())
        atoms_b = list(self.structure2.get_atoms())
        if len(atoms_a) > len(atoms_b):
            atoms_a = atoms_a[:len(atoms_b)]
        else:
            atoms_b = atoms_b[:len(atoms_a)]

        self.si.set_atoms(atoms_a, atoms_b)

        return self.si

    def getRMSD(self):

        superimpose = self.SuperimposeStructures()
        rmsd = superimpose.rms

        return rmsd


structures = {}
alphabet = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz"[::-1]
i = 0
for f in files_list:
    structures[f] = GetStructures(f)
    i += 2

superimposed_chains = Superimpose(list(structures.values())[0], list(structures.values())[1])
rmsd = superimposed_chains.getRMSD()
print(rmsd)
io = PDBIO()

