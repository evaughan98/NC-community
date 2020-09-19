import sys
import os
import json
from tabulate import tabulate

import RNA



# parameters ########################################

# command line arguments
L = int(sys.argv[1]) # sequence length
genotypes_section_total_number = int(sys.argv[2]) # total number of sections the genotype list is split
genotypes_section_number = int(sys.argv[3]) # section number considered in this run (starting from 1)

# fixed
alphabet = ["A","C","G","U"] # RNA alphabet



# directory and file structures ########################################

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

# help function (base transfer): genotype sequence -> genotype index
def gen_seq_to_gen_index(gen_seq,alphabet,L):
    gen_index = 0
    for site_index in range(L):
       gen_index += alphabet.index(gen_seq[site_index])*len(alphabet)**(L-site_index-1)
    return gen_index

# help function (base transfer): genotype index -> genotype sequence
def gen_index_to_gen_seq(gen_index,alphabet,L):
    gen_seq = []
    for site_index in range(L):
        value = gen_index // len(alphabet)**(L-site_index-1)
        letter = alphabet[value]
        gen_seq.append(letter)
        gen_index -= value*len(alphabet)**(L-site_index-1)
    return ''.join(gen_seq)

# analyse phenotypes for section of genotype list
phen_seqs = []
gp_map = {} # dictionary: genotype index -> phenotype index (with respect to the considered section of genotype list)

phen_seqs.append(''.join(["."]*L)) # set undefined phenotype to phenotype index 0

genotypes_section_length = n_genotypes / genotypes_section_total_number
genotypes_section_min_index = (genotypes_section_number-1)*genotypes_section_length
genotypes_section_max_index = (genotypes_section_number)*genotypes_section_length
if genotypes_section_number==genotypes_section_total_number:
    genotypes_section_max_index=n_genotypes # last section can be longer if n_genotypes / genotypes_section_total_number not an integer

for gen_index in xrange(genotypes_section_min_index,genotypes_section_max_index):
    gen_seq = gen_index_to_gen_seq(gen_index,alphabet,L)
    (structure, mfe) = RNA.fold(gen_seq)
    if structure not in phen_seqs:
        phen_seqs.append(structure)
    phen_index = phen_seqs.index(structure)
    if phen_index!=0: # ignore "undefined" phenotype (in order to save storage space)
        gp_map[gen_index] = phen_index



# write results ########################################

# write dictionary: genotype index -> phenotype index to .json file
file_name = results_directory_name + "L.%s_dictionary_genotype_index_to_phenotype_index_section.%s.json" % (L,genotypes_section_number)
with open(file_name,'wb') as handle:
    json.dump(gp_map,handle,sort_keys=True)

# write phenotypes to .txt file
table = [[i,phen_seqs[i]] for i,x in enumerate(phen_seqs)]
file_name = results_directory_name + "L.%s_phenotypes_section.%s.txt" % (L,genotypes_section_number)
f = open(file_name,"w")
f.write(tabulate(table,headers=["phen. index","phen. seq."]))
f.close()



print "done!"

quit()
