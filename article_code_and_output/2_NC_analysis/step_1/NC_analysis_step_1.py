import sys
import os
import numpy as np
import json
from tabulate import tabulate



# parameters ########################################

# command line arguments
L = int(sys.argv[1]) # sequence length
genotypes_section_total_number = int(sys.argv[2]) # total number of sections the genotype list is split (in "1_phenotype_analysis")

# fixed
alphabet = ["A","C","G","U"] # RNA alphabet



# directory and file structures ########################################

# input directory name
current_directory_name = os.getcwd()
current_directory_name_oneback = os.path.dirname(current_directory_name)
current_directory_name_twoback = os.path.dirname(current_directory_name_oneback)
input_directory_name = current_directory_name_twoback + "/1_phenotype_analysis/results/L.%s/" % (L)

# results directory name
results_directory_name = "results"
results_directory_name += "/L.%s" % (L)
results_directory_name += "/"

# create results directory (if not exist)
try:
    os.makedirs(results_directory_name)
except OSError:
    pass



# program ########################################

# total number of genotypes
n_genotypes = len(alphabet)**L

phen_seqs = []
phen_freqs = []
phen_maps = [] # (will contain:) dictionaries: genotype index -> phenotype index

# set undefined phenotype to phenotype index 0
phen_seqs.append(''.join(["."]*L))
phen_freqs.append(0)
phen_maps.append({})

# go through input (genotype list) sections:
for genotypes_section_number in range(1,genotypes_section_total_number+1):

    # load input
    input_file_name_1 = input_directory_name + "L.%s_dictionary_genotype_index_to_phenotype_index_section.%s.json" % (L,genotypes_section_number)
    input_file_name_2 = input_directory_name + "L.%s_phenotypes_section.%s.txt" % (L,genotypes_section_number)
    
    # dictionary: genotype index -> (initial) phenotype index
    with open(input_file_name_1, 'r') as handle:
        section_gp_map = json.load(handle)
    
    # (initial) phenotype indices, (initial) phenotype sequences
    section_phen_indices = list(np.loadtxt(input_file_name_2, usecols=(0,), dtype=int, skiprows=2, unpack=True))
    section_phen_seqs = list(np.loadtxt(input_file_name_2, usecols=(1,), dtype=str, skiprows=2, unpack=True))

    genotypes_section_length = n_genotypes / genotypes_section_total_number
    genotypes_section_min_index = (genotypes_section_number-1)*genotypes_section_length
    genotypes_section_max_index = (genotypes_section_number)*genotypes_section_length
    if genotypes_section_number==genotypes_section_total_number:
        genotypes_section_max_index=n_genotypes # last section can be longer if n_genotypes / genotypes_section_total_number not an integer

    for gen_index in xrange(genotypes_section_min_index,genotypes_section_max_index):
        if "%s" % (gen_index) in section_gp_map:
            phen_seq = section_phen_seqs[section_phen_indices.index(section_gp_map["%s" % (gen_index)])]
            if phen_seq not in phen_seqs:
                phen_seqs.append(phen_seq)
                phen_freqs.append(0)
                phen_maps.append({})
            phen_index = phen_seqs.index(phen_seq)
            phen_freqs[phen_index] += 1
            phen_maps[phen_index][gen_index] = phen_index
        else: # undefined phenotype
            phen_freqs[0] += 1



# write results ########################################

# write phenotypes to .txt file
table = [[i,phen_seqs[i],phen_freqs[i]] for i,x in enumerate(phen_seqs)]
file_name = results_directory_name + "L.%s_phenotypes.txt" % (L)
f = open(file_name,"w")
f.write(tabulate(table,headers=["phen. index","phen. seq.","phen. freq."]))
f.close()

# write dictionaries: genotype index -> phen index to .json file
for phen_index,phen_map in enumerate(phen_maps):
    if phen_index!=0: # exclude undefined phenotype
        file_name = results_directory_name + "L.%s_dictionary_genotype_index_to_phenotype_index.%s.json" % (L,phen_index)
        with open(file_name,'wb') as handle:
            json.dump(phen_map,handle,sort_keys=True)



print "done!"

quit()
