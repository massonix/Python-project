# all functions and classes 1


# multifunctions.py contains all the functions and classes required by multicomplex.py to
# reconstruct a complete complex from a set of interacting pairs.

from Bio.PDB import *
from Bio import pairwise2
import os


def GetStructures(pdbfile, ids):
    """
    Returns a structure object given a pdb file.

    Args:
        pdbfile: a string describing the path of a pdb file.
        ids: a string that will be used to assign new chain ids
            to the structure object.

    Returns:
        A structure object from the Bio.PDB package, which follows
        the SMCRA (Structure/Model/Chain/Residue/Atom) architecture.

    """

    parser_pdb = PDBParser()
    structure = parser_pdb.get_structure(pdbfile[0:-4], pdbfile)
    current_ids = [x.id for x in structure[0].get_chains()]
    i = 0

    for chain in structure[0].get_chains():
        if current_ids[0] == ids[1] or current_ids[1] == ids[0]:
            list(structure[0].get_chains())[0].id = 1000
            list(structure[0].get_chains())[1].id = 1001
            list(structure[0].get_chains())[0].id = ids[0]
            list(structure[0].get_chains())[1].id = ids[1]
            break
        elif ids[i] != current_ids[i]:
            chain.id = ids[i]
        i += 1
    return structure


def get_interactions(list_atoms1, list_atoms2, dist):
    """
    Get all the pairs of interacting residues from two proteins that
    are within a certain distance.

    Args:
        list_atoms1: list of all atoms from the first protein.
        list_atoms2: list of all atoms from the second protein. Each
                     atom should be an atom object from the Bio.PDB
                     package, following the SMCRA architecture.
        dist: maximum distance (in Armstrongs) at which two beta carbons
              from the two proteins need to be in order to be returned.

    Returns:
        a list of tuples, each of which contains two residue ids corresponding
        to two residues from both proteins that are within dist Armstrongs.

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
    Returns a sequence alignment given 2 structures.

    Args:
        structure1: first structure, model or chain object from
                    Bio.PDB, following the SMCRA architecture.
        structure2: second structure, model or chain object from
                    Bio.PDB, following the SMCRA architecture.

    Returns:
        A list containing the following elements:
        1. and 2. Strings showing the alignment alignment between
                 both structures.
        3. Score of the alignment.
        4. and 5. Begin and end (where the alignment occurs).

    """

    ppb = PPBuilder()
    for pp in ppb.build_peptides(structure1):
        sequence1 = pp.get_sequence()
    for pp in ppb.build_peptides(structure2):
        sequence2 = pp.get_sequence()

    alignment = pairwise2.align.globalxx(sequence1, sequence2)
    return alignment


class Superimpose(object):
    """
    Superimpose objects allow the superimposition of two structures.

    Arguments:
        structure1: first structure, model or chain object from
                    Bio.PDB, following the SMCRA architecture.
        structure2: second structure, model or chain object from
                    Bio.PDB, following the SMCRA architecture.

    """
    def __init__(self, structure1, structure2):
        """Inits Superimpose class with structure1 and structure2"""

        self.structure1 = structure1
        self.structure2 = structure2
        self.si = Superimposer()

    def SuperimposeStructures(self):
        """
        Returns the structure object resulting from the structural superimposition
        of the two structures passed as attributes when defining a new instance of
        Superimpose.
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
        """
        Returns the rmsd of the structural superimposition between structure1 and
        structure2.
        """

        superimpose = self.SuperimposeStructures()
        rmsd = superimpose.rms

        return rmsd


models = []
models2 = []
j = 0

def impose_clash(str1, strs, k, i, num_contact, c, similar_chains, stoichiometry, outdir, max_files):
    """
    impose_clash is the core recursive function of our project. It takes a structure as input (str1) and tries to add
    a new chain out of several structures (strs) to build a larger complex. Each new potential complex needs to meet
    the following criteria in order to call again the function with the resulting complex as input:
    1. Similarity: As the superimposition is based on sequence similarity, both structures need to have two similar chains
        in common in order for the function to superimpose them.
    2. Clashes: once the superimposition is made, the structure is feasible only if one structure does not clash with the
       other (i.e. the alpha carbons are not in close contact between chains).
    3. Duplicities: as the function is recursive, some equal complexes will be achieved by different paths. The functions
        calls itself again if and only if the current complex has not been built before.
    4. Stoichiometry: if the desired stoichiometry is met (i.e. A2B2), the function will save the complex in a file and
        return None; if the current complex exceeds the desired stoichiometry (i.e. A3B1), the function will return None.
    5. Recursion depth: computed as the sum of all the chains in the stoichiometry (i.e. k = 4 if stoich = A2B2; k=7 if
       stoich = A7). If the current complex has a # of chains equal to k without meeting the stoichiometry, the function
       will return None.
    6. Number of output files: for some purposes the script may take too long to complete all recursions. In such
       scenarios, one can limit the number of desired complexes he or she wants to obtain.


    Args:
        str1: the structure object that the function will try to superimpose against all other structures.
        strs:
        k: recursion depth
        i:
        num_contact: number of contacts allowed between chains, considering a contact two beta carbons from different
                     chains closer than 2A.
        c:
        similar_chains: dictionary that has as keys all the chain identfiers from all structures in strs, and as values
                        a chain id to which they are similar (all similar chains map to the same chain id).
        stoichiometry: a dictionary describing the desired stoiciometry. For example, {"A":3, "B":1} if the stoichiometry
                        passed by the user is "A3B1".

    Returns:
        The function itself returns None, but modifies the global variables models1 and saves the output complexes as
        pdb files.
    """

    global models
    global models2
    global j

    alphabet = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz"

    if max_files:
        current_file_num = len(list(filter(lambda x: x.startswith("complex"), os.listdir(outdir))))
        if current_file_num >= max_files:
            return

    current_stoich = {}
    stoich_list = [similar_chains[x.id] for x in str1.get_chains()]
    stoich_set = set(stoich_list)
    stoich_list2 = [stoich_list.count(x) for x in stoich_set]
    stoich_list2.sort(reverse = True)

    for index, s in enumerate(stoich_list2):
        current_stoich[alphabet[index]] = s
        if current_stoich[alphabet[index]] > stoichiometry[alphabet[index]]:
            return

    if  current_stoich == stoichiometry:
        for ind, ch in enumerate(list(str1.get_chains())):
            ch.id = alphabet[ind]
        models2.append(str1)
        print(models2)
        j += 1
        io = PDBIO()
        io.set_structure(str1.copy()[0])
        io.save(outdir + "/complex" + str(j) + ".pdb")
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
                    clashes = get_interactions(list(chain_diff2.get_atoms()), list(str3[0].get_atoms()), 2)

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
                            models.append(str3)
                            impose_clash(str3, strs, k, i + 1, num_contact, c, similar_chains, stoichiometry, outdir, max_files)
                            if max_files:
                                current_file_num = len(
                                    list(filter(lambda x: x.startswith("complex"), os.listdir(outdir))))
                                if current_file_num >= max_files:
                                    return
                            elif i >= k:
                                return

                else:

                    fails += 1

    if fails == len(strs) * len(chains1):
        return
