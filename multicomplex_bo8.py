##stoichiometry + indefinite recurssion + duplicities removal

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


def GetStructures(pdbfile, ids):
    """
    Given a pdbfile, gets it's structure.
    :param pdbfile: pdb file to open
    :return: structure object
    """
    parser = PDBParser()
    structure = parser.get_structure(pdbfile[0:-4], pdbfile)
    current_ids = [x.id for x in structure[0].get_chains()]
    i = 0

    for chain in structure[0].get_chains():
        if ids[i] != current_ids[i]:
            chain.id = ids[i]
        i += 1
    return structure


def get_interactions(list_atoms1, list_atoms2, dist):
    """given 2 lists of atoms corresponding to 2 different chains,
    returns a tuple with 3 elements:
        1. tuple of chains that interact, i.e. ("A","B")
        2. tuple of residue number that interact, i.e. (125, 543)
        3. distance
    """
    beta_carbons1 = list(filter(lambda x: x.get_id() == "CB", list_atoms1))
    beta_carbons2 = list(filter(lambda x: x.get_id() == "CB", list_atoms2))
    ns = NeighborSearch(beta_carbons1)
    interactions = []

    for atom in beta_carbons2:
        interact = ns.search(atom.get_coord(), dist)
        interactions.extend(
            [tuple(sorted([str(atom.get_parent().resname), str(x.get_parent().resname)])) for x in interact])
    return interactions


def Alignsequence(structure1, structure2):
    """
    Given 2 structures, gets the sequence and aligns it
    :param structure1: structure 1 to align
    :param structure2: structure 2 to align
    :return: Alignment
    """
    ppb = PPBuilder()
    for pp in ppb.build_peptides(structure1):
        sequence1 = pp.get_sequence()
    for pp in ppb.build_peptides(structure2):
        sequence2 = pp.get_sequence()

    alignment = pairwise2.align.globalxx(sequence1, sequence2)
    return alignment

models = []
models2 = []
def impose_clash(str1, strs, k, i, num_contact, c, similar_chains, stoichiometry):
    """given 2 structures, impose_clash rotates the second structure to superimpose with the common chain. If a
    clash if found (aa from the other chain closer than 2 amstrongs), structure1 is returned. If no clash is found,
    the superimposition is returned
    """
    global models
    global models2
    alphabet = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz"

    current_stoich = {}
    print(stoich_dict)
    stoich_list = [similar_chains[x.id] for x in str1.get_chains()]
    stoich_set = set(stoich_list)
    stoich_list2 = [stoich_list.count(x) for x in stoich_set]
    stoich_list2.sort(reverse = True)
    for index, s in enumerate(stoich_list2):
        current_stoich[alphabet[index]] = s
        if current_stoich[alphabet[index]] > stoich_dict[alphabet[index]]:
            return
    if  current_stoich == stoich_dict:
        for ind, ch in enumerate(list(str1.get_chains())):
            ch.id = alphabet[ind]
        models2.append(str1)
        io = PDBIO()
        io.set_structure(str1)
        io.save("selected_pdbs3/dimer" + str(id(str1)) + ".pdb")
        return

    elif i >= k:
        return

    fails = 0
    chains1 = list(str1.get_chains())


    for str2 in strs:
        chains2 = list(str2.get_chains())
        for chain1 in chains1:
            for chain2 in chains2:
                str3 = str1.copy()
                str4 = str2.copy()
                copies3 = dict([(x, y) for x, y in zip(chains1, list(str3[0].get_chains()))])
                copies4 = dict([(x, y) for x, y in zip(chains2, list(str4[0].get_chains()))])


                if similar_chains[chain1.id] == similar_chains[chain2.id]:
                    common_chain1 = copies3[chain1]
                    common_chain2 = copies4[chain2]
                    superimposed_chains = Superimpose(common_chain1, common_chain2)
                    superimposed_chains_fin = superimposed_chains.SuperimposeStructures()
                    superimposed_chains_fin.apply(list(str4[0].get_atoms()))
                    c += 1
                    chain_diff = [x for x in str4[0].get_chains() if x.id != common_chain2.id]
                    chain_diff2 = chain_diff[0].copy()
                    chain_diff2.id = id(chain_diff2) + c
                    clashes = get_interactions(list(chain_diff2.get_atoms()), list(str3[0].get_atoms()), 1.5)
                    if len(clashes) >= num_contact:
                        fails += 1

                    else:
                        str3[0].add(chain_diff2)
                        similar_chains[chain_diff2.id] = similar_chains[str2[0][chain_diff[0].get_id()].id]
                        repeated = False
                        str5 = str3.copy()
                        for model in models:
                            superimposed_models = Superimpose(model, str5[0])
                            rmsd = superimposed_models.getRMSD()
                            if rmsd < 10 and len(list(model.get_chains())) == len(list(str5.get_chains())):
                                repeated = True
                        if not repeated:
                            str5.id = str(id(str1))
                            models.append(str3)
                            impose_clash(str3, strs, k, i + 1, 10, c, similar_chains, stoich_dict)


                else:

                    fails += 1


    if fails == len(strs) * len(chains1):
        return


structures = {}
alphabet = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz"[::-1]
i = 0
for f in files_list:
    structures[f] = GetStructures(f, alphabet[i:i+2])
    i += 2

stoich_in = "A4"
k = sum([int(x) for x in list(stoich_in[1::2])])
stoich_dict = {stoich_in[i]: int(stoich_in[i+1]) for i in range(0, len(stoich_in) - 1, 2)}

"""get the similar chains and remove those dimers that do not interact with anyone else
The structures dict contains the filenames as keys and the structure objects as values.
After the following chunk of code we remove those files:str that do not have a similar chain with
another dimer."""

structures2 = structures.copy()
scores = {}
all_files = []
discarded_files = set()
superimpositions = []
structures_iter = list(structures2.items())
similar_chains = {}
i = 0
for index1, items1 in enumerate(structures_iter):
    file1 = items1[0]
    structure1 = items1[1]
    chains1 = list(structure1[0].get_chains())
    for file2, structure2 in structures_iter[index1:len(structures_iter)]:
        chains2 = list(structure2[0].get_chains())

        for chain1 in chains1:
            i += 1
            for chain2 in chains2:
                i += 1
                Alignment = Alignsequence(chain1, chain2)
                score = Alignment[0][2] / len(Alignment[0][0])
                if score > 0.95:
                    similar_chains.setdefault(chain2.id, chain1.id)

    if chains2[0].id in similar_chains and chains2[1].id in similar_chains:
        if similar_chains[chains2[0].id] == similar_chains[chains2[0].id] and similar_chains[chains2[1].id] == similar_chains[
            chains2[1].id]:
            structures.pop(file1)
            del similar_chains[chains2[0].id]
            del similar_chains[chains2[1].id]

strs = list(structures.values())
print(strs)
impose_clash(strs[0], strs, k, 2, 10, 0, similar_chains, stoich_dict)

#print(list(models[-1].get_chains()))
#for model in models:
 #   io = PDBIO()
#    io.set_structure(model)
 #   io.save("selected_pdbs3/dimer_prova" + str(i) +".pdb")
  #  i += 1
# # i = 0
# # io = PDBIO()


# # for i in range(len(strs_out)):
# #     with open("./selected_pdbs2/dimer" + str(i+1) + ".pdb", "w") as fp:

# #         io.set_structure(strs_out[i][0])
# #         io.save(fp)




