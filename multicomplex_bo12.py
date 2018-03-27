# With argparse, max_number of files working, module stratification, documented.

from Bio.PDB import *
from Bio import pairwise2
from multifunctions import *
import sys
import os
import argparse
import string


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Reconstructs the complete macrocomplex given a set of interacting pairs (prot-prot).")

    parser.add_argument('-i', '--input',
                        dest="infiles",
                        action="store",
                        default=None,
                        help="Input PDB files or directory with PDB files.")

    parser.add_argument('-o', '--output',
                        dest="outdir",
                        action="store",
                        default=None,
                        help="Directory where the macrocomplexes will be saved as pdbs.")

    parser.add_argument('-v', '--verbose',
                        dest="verbose",
                        action="store_true",
                        default="False",
                        help="Print log in stderr")

    parser.add_argument('-st', '--stoich',
                        dest="stoichiometry",
                        action="store",
                        default=None,
                        type=str,
                        help="String with the stoichiometry of the Macrocomplex.")

    parser.add_argument('-opt', '--optimization',
                        dest="optimization",
                        action="store_true",
                        default="False",
                        help="Optimize the output")

    parser.add_argument('-cn', '--contact_num',
                        dest="contact_num",
                        action="store",
                        default=None,
                        type=int,
                        help="Maximum number of permited contacts")

    parser.add_argument('-f', '--files_num',
                        dest="max_files",
                        action="store",
                        default=None,
                        type=int,
                        help="Number of files to output")

    options = parser.parse_args()

    # 1. Save all pdb files passed by the user in a list: files_list.

    inter = options.infiles
    if os.path.isdir(inter):
        files_list = []
        for files in os.listdir(inter):
            if files.endswith(".pdb"):
                files_list.append(os.path.join(inter, files))



    # 2. Create a structures dictionary, with the files path as key and the structure object as values.

    structures = {}
    alphabet = string.printable
    i = 0
    for f in files_list:
        structures[f] = GetStructures(f, alphabet[i:i+2])
        i += 2


    # 3. Process the stoichiometry input to calculate the recurssion depth (k) and create the dictionary stoich_dict.

    stoich_in = options.stoichiometry
    stoich_dict = {}
    i = 0
    j = 0
    while i < len(stoich_in):
        if stoich_in[i].isalpha():
            stoich_dict[stoich_in[i]] = ""
            j += 1
            while stoich_in[j].isdigit():
                stoich_dict[stoich_in[i]] += str(stoich_in[j])
                if j < (len(stoich_in)-1):
                    j += 1
                else:
                    break
        i += 1
        j = i

    stoich_dict = {x:int(y) for x,y in stoich_dict.items()}
    k = sum(list(stoich_dict.values()))


    # 4. Map all the similar chains (>95% sequence similarity): similar_chains

    structures2 = structures.copy()
    scores = {}

    structures_iter = list(structures2.items())
    similar_chains = {}
    i = 0
    print(len(structures))
    for index1, items1 in enumerate(structures_iter):
        file1 = items1[0]
        structure1 = items1[1]
        chains1 = list(structure1[0].get_chains())
        for file2, structure2 in structures_iter[index1:len(structures_iter)]:
            chains2 = list(structure2[0].get_chains())

            for chain1 in chains1:
                i += 1
                for chain2 in chains2:
                    if chain2.id in similar_chains:
                        continue
                    i += 1
                    Alignment = Alignsequence(chain1, chain2)
                    score = Alignment[0][2] / len(Alignment[0][0])
                    if score > 0.95:
                        similar_chains.setdefault(chain2.id, chain1.id)


    # 5. Remove those structures that do not share any similar chain with another structure and thus cannot be superimposed.

    similar_chains_keys = sorted(list(similar_chains.keys()), reverse = True)
    similar_chains_values = [similar_chains[x] for x in similar_chains_keys]
    similar_chains_values_iter = [(similar_chains_values[x],similar_chains_values[x+1]) for x in range(0, len(similar_chains_values), 2)]

    for i, chs in enumerate(similar_chains_values_iter):
        if similar_chains_values.count(chs[0]) == 1 and similar_chains_values.count(chs[1]) == 1:
            structures.pop(files_list[i])
        elif chs[0] == chs[1] and similar_chains_values.count(chs[0]) == 2:
            structures.pop(files_list[i])

    del structures2
    del scores
    del structures_iter
    del similar_chains_keys
    del similar_chains_values
    del similar_chains_values_iter


    # 6. Call the recursive function and built the models with the desired stoichiometry.

    if options.contact_num:
        contacts = options.contact_num
    else:
        contacts = 5

    strs = list(structures.values())
    impose_clash(strs[0], strs, k, 2, contacts, 0, similar_chains, stoich_dict, options.outdir, options.max_files)


