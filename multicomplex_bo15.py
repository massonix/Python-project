# Works with proteasome, phosphate dehydratase and nucleosome. Only thing left: ensure removes duplications, finish documenting
# make it installable, write README.

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
    alphabet = list(range(1000))
    i = 0
    for f in files_list:
        structures[f] = GetStructures(f, alphabet[i:i+2])
        i += 2

    if options.verbose:
        sys.stderr.write("%d PDB files found: \n" % len(structures))
        for f in structures.keys():
            sys.stderr.write(f + "\n")

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

    if options.verbose:
        sys.stderr.write("Stoichiometry: %s\n" % options.stoichiometry)


    # 4. Map all the similar chains (>95% sequence similarity): similar_chains

    structures2 = structures.copy()
    scores = {}

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
                    if chain2.id in similar_chains:
                        continue
                    i += 1
                    Alignment = Alignsequence(chain1, chain2)
                    #print([x.resname for x in chain1.get_residues()])
                    #print([x.resname for x in chain2.get_residues()])
                    score = Alignment[0][2] / len(Alignment[0][0])
                    if score > 0.95:
                        similar_chains.setdefault(chain2.id, chain1.id)


    # 5. Remove those structures that do not share any similar chain with another structure and thus cannot be superimposed.

    similar_chains_keys = sorted(list(similar_chains.keys()))
    similar_chains_values = [similar_chains[x] for x in similar_chains_keys]
    similar_chains_keys_iter = [(similar_chains_keys[x],similar_chains_keys[x+1]) for x in range(0, len(similar_chains_keys), 2)]
    similar_chains_values_iter = [(similar_chains_values[x],similar_chains_values[x+1]) for x in range(0, len(similar_chains_values), 2)]

    for i, chs in enumerate(similar_chains_values_iter):
        if similar_chains_values.count(chs[0]) == 1 and similar_chains_values.count(chs[1]) == 1:
            structures.pop(files_list[i])
            del similar_chains[similar_chains_keys_iter[i][0]]
            del similar_chains[similar_chains_keys_iter[i][1]]
        elif chs[0] == chs[1] and similar_chains_values.count(chs[0]) == 2:
            structures.pop(files_list[i])
            del similar_chains[similar_chains_keys_iter[i][0]]
            del similar_chains[similar_chains_keys_iter[i][1]]

    if not similar_chains:
        raise ValueError("Unable to superimpose: no common chains among structures")

    del structures2
    del scores
    del structures_iter
    del similar_chains_keys
    del similar_chains_values
    del similar_chains_keys_iter


    # 6. Check if there are enough different chains to achieve the desired stoichiometry.

    stoich_set = set(list(similar_chains.values()))

    if len(stoich_set) < len(stoich_dict):
        raise ValueError("Impossibe stoichiometry: The provided stoichiometry contains %d different chains, but "
                         "the input PDBs only have %d unique chains." % (len(stoich_dict), len(stoich_set)))


    # 7. Call the recursive function and build the models with the desired stoichiometry.

    if options.contact_num:
        contacts = options.contact_num
    else:
        contacts = 5

    strs = [structures[x] for x in files_list]
    stoich_list = list(stoich_set)
    strs_similar_dic = {strs[x]:similar_chains_values_iter[x] for x in range(len(strs))}


    for ind, uniq_chain1 in enumerate(stoich_list):
        for uniq_chain2 in stoich_list[ind:]:
            dimer = [x for x in strs if strs_similar_dic[x] == (uniq_chain1, uniq_chain2) or strs_similar_dic[x] == (uniq_chain2, uniq_chain1)]
            if dimer:
                if uniq_chain1 == uniq_chain2:
                    if len(stoich_dict) == 1:
                        strs2 = [x for x in strs if strs_similar_dic[x] == (uniq_chain1, uniq_chain2)]
                        if options.verbose:
                            sys.stderr.write("Starting recursion...\n")
                        impose_clash(dimer[0], strs2, k, 2, contacts, 0, similar_chains, stoich_dict, options.outdir,
                        options.max_files, options.verbose)
                    else:
                        continue
                else:
                    if len(stoich_dict) == 1:
                        continue
                    else:
                        strs2 = [x for x in strs if (uniq_chain1 in strs_similar_dic[x] or uniq_chain2 in strs_similar_dic[x]) and
                                 strs_similar_dic[x] != (uniq_chain1, uniq_chain1) and strs_similar_dic[x] != (uniq_chain2, uniq_chain2)]
                        if options.verbose:
                            sys.stderr.write("Starting recursion...\n")
                        impose_clash(dimer[0], strs, k, 2, contacts, 0, similar_chains, stoich_dict, options.outdir,
                                     options.max_files, options.verbose)


    if options.verbose:
        sys.stderr.write("Program finished correctly. \n")




