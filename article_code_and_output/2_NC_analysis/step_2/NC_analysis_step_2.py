import sys
import os
import numpy as np
import json
from tabulate import tabulate



# parameters ########################################

# command line arguments
L = int(sys.argv[1]) # sequence length
phen_index = int(sys.argv[2]) # phenotype index

# fixed
alphabet = ["A","C","G","U"] # RNA alphabet



# directory and file structures ########################################

# input directory name
current_directory_name = os.getcwd()
current_directory_name_oneback = os.path.dirname(current_directory_name)
input_directory_name = current_directory_name_oneback + "/step_1/results/L.%s/" % (L)

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

# load phenotype map: dictionary: genotype index -> phenotype index
input_file_name = input_directory_name + "L.%s_dictionary_genotype_index_to_phenotype_index.%s.json" % (L,phen_index)

with open(input_file_name, 'r') as handle:
    phen_map = json.load(handle)



## find NCs ########################################

# help function: find neutral neighbour genotype indices for a genotype index
def find_neutral_nb_gen_indices(gen_index,phen_map,alphabet,L):
    gen_seq = gen_index_to_gen_seq(gen_index,alphabet,L)
    gen_seq_array = list(gen_seq)
    neutral_nb_gen_indices = set([]) # set of all neutral neighbour genotype indices
    # go through all characters of the genotype sequence
    for i,character in enumerate(gen_seq_array):
        # go through all letters of the alphabet
        for letter in alphabet:
            # mutate character to "value" of letter if character != letter
            if character!=letter:
                mut_gen_seq_array = [j for j in gen_seq_array] # to preserve "gen_seq_array"
                # mutate
                mut_gen_seq_array[i] = letter
                mut_gen_index = gen_seq_to_gen_index(''.join(mut_gen_seq_array),alphabet,L)
                # if neutral mutation, add mutated genotype index to set of neutral neighbour genotype indices
                if "%s" % (mut_gen_index) in phen_map:
                    neutral_nb_gen_indices.add(mut_gen_index)
    return neutral_nb_gen_indices

# find NCs (adapted from S. Schaper, PhD thesis, University of Oxford, 2012)
NCs = []
tested_gen_indices = set([])

for gen_index in phen_map:
    gen_index = int(gen_index)
    # if not already part of a NC, determine the "new" NC of this genotype
    if gen_index not in tested_gen_indices:
        U = set([gen_index]) # set of genotype indices that are part of NC, but for which one-point mutational neighbourhood is not yet checked
        V = set([]) # set of genotype indices that are part of NC, and for which one-point mutational neighbourhood has been checked
        while len(U) > 0:
            gen1_index = U.pop()
            neutral_nb_gen_indices = find_neutral_nb_gen_indices(gen1_index,phen_map,alphabet,L)
            for gen2_index in neutral_nb_gen_indices:
                if gen2_index not in V:
                    U.add(gen2_index)
            V.add(gen1_index)
        NC = sorted(list(V))
        NCs.append(NC)
        for tested_gen_index in NC:
            tested_gen_indices.add(tested_gen_index)



## analyse NC characteristics ########################################

# help function: find neutral mutations per site for a genotype index
def find_neutral_mut_per_site(gen_index,phen_map,alphabet,L):
    gen_seq = gen_index_to_gen_seq(gen_index,alphabet,L)
    gen_seq_array = list(gen_seq)
    gen_neutral_mut_per_site = [0]*L
    gen_acc_alt_phen_indices = set([])
    # go through all characters of the genotype sequence
    for i,character in enumerate(gen_seq_array):
        # go through all letters of the alphabet
        for letter in alphabet:
            # mutate character to "value" of letter if character != letter
            if character!=letter:
                mut_gen_seq_array = [j for j in gen_seq_array] # to preserve "gen_seq_array"
                # mutate
                mut_gen_seq_array[i] = letter
                mut_gen_index = gen_seq_to_gen_index(''.join(mut_gen_seq_array),alphabet,L)
                # if neutral mutation, add it to the specific site
                if "%s" % (mut_gen_index) in phen_map:
                    gen_neutral_mut_per_site[i] += 1
    return gen_neutral_mut_per_site

# analyse NCs
gNC_map = {} # dictionary: genotype index -> NC index

NC_avg_neutral_mut_per_site_list = []
NC_SD_neutral_mut_per_site_list = []

NC_sizes = []
NC_robs = []

for NC_index,NC in enumerate(NCs):

    gen_neutral_mut_per_site_storage = []

    for gen_index in NC:
    
        gNC_map[gen_index] = NC_index # fill dictionary
    
        gen_neutral_mut_per_site = find_neutral_mut_per_site(gen_index,phen_map,alphabet,L)
    
        gen_neutral_mut_per_site_storage.append(gen_neutral_mut_per_site)

    # NC neutral mutations per site analysis
    NC_avg_neutral_mut_per_site = [0 for l in range(L)]
    NC_SD_neutral_mut_per_site = [0 for l in range(L)]

    for k in range(len(NC)):
        for l in range(L):
            NC_avg_neutral_mut_per_site[l] += gen_neutral_mut_per_site_storage[k][l]
    for l in range(L):
        NC_avg_neutral_mut_per_site[l] = float(NC_avg_neutral_mut_per_site[l])/len(NC)

    for k in range(len(NC)):
        for l in range(L):
            NC_SD_neutral_mut_per_site[l] += (gen_neutral_mut_per_site_storage[k][l]-NC_avg_neutral_mut_per_site[l])**2
    for l in range(L):
        NC_SD_neutral_mut_per_site[l] = np.sqrt(float(NC_SD_neutral_mut_per_site[l])/len(NC))
    
    NC_avg_neutral_mut_per_site_list.append(NC_avg_neutral_mut_per_site)
    NC_SD_neutral_mut_per_site_list.append(NC_SD_neutral_mut_per_site)

    # NC characteristics
    NC_size = len(NC)
    NC_rob = float(sum(NC_avg_neutral_mut_per_site))/L/(len(alphabet)-1)

    NC_sizes.append(NC_size)
    NC_robs.append(NC_rob)



# write results ########################################

# write dictionary: genotype index -> NC index to .json file
file_name = results_directory_name + "L.%s_dictionary_genotype_index_to_NC_index_phenotype_index.%s.json" % (L,phen_index)
with open(file_name,'wb') as handle:
    json.dump(gNC_map,handle,sort_keys=True)

# write NC characteristics to .txt file
table = [[i,NC_sizes[i],NC_robs[i]] for i,x in enumerate(NCs)]
file_name = results_directory_name + "L.%s_neutral_component_characteristics_phenotype_index.%s.txt" % (L,phen_index)
f = open(file_name,"w")
f.write(tabulate(table,headers=["NC index","NC size","NC rob."]))
f.close()

# write NC average neutral mutations per site to .txt file
table = [[i] for i,x in enumerate(NCs)]
for i,x in enumerate(NCs):
    for l in range(L):
        table[i].append(NC_avg_neutral_mut_per_site_list[i][l])
file_name = results_directory_name + "L.%s_neutral_component_average_neutral_mutations_per_site_phenotype_index.%s.txt" % (L,phen_index)
f = open(file_name,"w")
headers=["NC index"] + [""]*L
f.write(tabulate(table,headers=headers))
f.close()

# write NC SD neutral mutations per site to .txt file
table = [[i] for i,x in enumerate(NCs)]
for i,x in enumerate(NCs):
    for l in range(L):
        table[i].append(NC_SD_neutral_mut_per_site_list[i][l])
file_name = results_directory_name + "L.%s_neutral_component_SD_neutral_mutations_per_site_phenotype_index.%s.txt" % (L,phen_index)
f = open(file_name,"w")
headers=["NC index"] + [""]*L
f.write(tabulate(table,headers=headers))
f.close()



print "done!"

quit()
