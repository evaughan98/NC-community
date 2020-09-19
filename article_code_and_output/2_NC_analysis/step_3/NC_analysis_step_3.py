import sys
import os
import numpy as np
import json
from tabulate import tabulate



# parameters ########################################

# command line arguments
L = int(sys.argv[1]) # sequence length



# directory and file structures ########################################

# input directory names
current_directory_name = os.getcwd()
current_directory_name_oneback = os.path.dirname(current_directory_name)
input_directory_name_1 = current_directory_name_oneback + "/step_1/results/L.%s/" % (L)
input_directory_name_2 = current_directory_name_oneback + "/step_2/results/L.%s/" % (L)

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

## load input ########################################

# phenotype characteristics
input_file_name_1 = input_directory_name_1 + "L.%s_phenotypes.txt" % (L)

phen_indices = list(np.loadtxt(input_file_name_1, usecols=(0,), dtype=int, skiprows=2, unpack=True))
phen_seqs = list(np.loadtxt(input_file_name_1, usecols=(1,), dtype=str, skiprows=2, unpack=True))
phen_freqs = list(np.loadtxt(input_file_name_1, usecols=(2,), dtype=int, skiprows=2, unpack=True))

# NC characteristics and NC maps
NC_indices = [] # new NC indices
NC_phen_indices = []
NC_phen_seqs = []
NC_sizes = []
NC_robs = []
NC_avg_neutral_mut_per_site_list = []
NC_SD_neutral_mut_per_site_list = []

gNC_map = {} # dictionary: genotype index -> new NC index

for phen_index in phen_indices:
    if phen_index!=0: # exclude undefined phenotype
    
        input_file_name_2 = input_directory_name_2 + "L.%s_neutral_component_characteristics_phenotype_index.%s.txt" % (L,phen_index)
        try:
            load_NC_indices = list(np.loadtxt(input_file_name_2, usecols=(0,), dtype=int, skiprows=2, unpack=True)) # old NC indices
            load_NC_sizes = list(np.loadtxt(input_file_name_2, usecols=(1,), dtype=int, skiprows=2, unpack=True))
            load_NC_robs = list(np.loadtxt(input_file_name_2, usecols=(2,), skiprows=2, unpack=True))
        except TypeError: # only one NC
            load_NC_indices = [int(np.loadtxt(input_file_name_2, usecols=(0,), dtype=int, skiprows=2, unpack=True))]
            load_NC_sizes = [int(np.loadtxt(input_file_name_2, usecols=(1,), dtype=int, skiprows=2, unpack=True))]
            load_NC_robs = [float(np.loadtxt(input_file_name_2, usecols=(2,), skiprows=2, unpack=True))]
        
        input_file_name_3 = input_directory_name_2 + "L.%s_neutral_component_average_neutral_mutations_per_site_phenotype_index.%s.txt" % (L,phen_index)
        load_NC_average_neutral_mut_per_site_list = [[] for load_NC_index in load_NC_indices]
        for l in range(L):
            try:
                load = list(np.loadtxt(input_file_name_3, usecols=(1+l,), skiprows=2, unpack=True))
            except TypeError: # only one NC
                load = [np.loadtxt(input_file_name_3, usecols=(1+l,), skiprows=2, unpack=True)]
            for i,load_NC_index in enumerate(load_NC_indices):
                load_NC_average_neutral_mut_per_site_list[i].append(load[i])
        
        input_file_name_4 = input_directory_name_2 + "L.%s_neutral_component_SD_neutral_mutations_per_site_phenotype_index.%s.txt" % (L,phen_index)
        load_NC_SD_neutral_mut_per_site_list = [[] for load_NC_index in load_NC_indices]
        for l in range(L):
            try:
                load = list(np.loadtxt(input_file_name_4, usecols=(1+l,), skiprows=2, unpack=True))
            except TypeError: # only one NC
                load = [np.loadtxt(input_file_name_4, usecols=(1+l,), skiprows=2, unpack=True)]
            for i,load_NC_index in enumerate(load_NC_indices):
                load_NC_SD_neutral_mut_per_site_list[i].append(load[i])
        
        input_file_name_5 = input_directory_name_2 + "L.%s_dictionary_genotype_index_to_NC_index_phenotype_index.%s.json" % (L,phen_index)
        with open(input_file_name_5, 'r') as handle:
            load_gNC_map = json.load(handle) # dictionary: genotype index -> old NC index
            
        old_NC_index_new_NC_index_map = {} # dictionary: old NC index -> new NC index
        
        base = len(NC_indices)
        for old_NC_index in load_NC_indices:
        
            old_NC_index_new_NC_index_map[old_NC_index] = base+old_NC_index
            
            NC_indices.append(base+old_NC_index)
            NC_phen_indices.append(phen_index)
            NC_phen_seqs.append(phen_seqs[phen_index])
            NC_sizes.append(load_NC_sizes[old_NC_index])
            NC_robs.append(load_NC_robs[old_NC_index])
            NC_avg_neutral_mut_per_site_list.append(load_NC_average_neutral_mut_per_site_list[old_NC_index])
            NC_SD_neutral_mut_per_site_list.append(load_NC_SD_neutral_mut_per_site_list[old_NC_index])
            
        for gen_index in load_gNC_map:
            gen_index = int(gen_index)
            old_NC_index = load_gNC_map["%s" % (gen_index)]
            new_NC_index = old_NC_index_new_NC_index_map[old_NC_index]
            gNC_map[gen_index] = new_NC_index
            
        
        
## sort NC and phenotype characteristics and add ranks ########################################

# sort NC characteristics according to NC size
NC_sizes, NC_indices, NC_phen_indices, NC_phen_seqs, NC_robs, NC_avg_neutral_mut_per_site_list, NC_SD_neutral_mut_per_site_list = (list(x) for x in zip(*sorted(zip(NC_sizes,NC_indices,NC_phen_indices,NC_phen_seqs,NC_robs,NC_avg_neutral_mut_per_site_list,NC_SD_neutral_mut_per_site_list),reverse=True)))
NC_ranks = [i+1 for i,x in enumerate(NC_indices)] # rank, largest NC has rank 1

# sort phenotype characteristics according to phenotype frequency
phen_freqs[1:], phen_indices[1:], phen_seqs[1:] = (list(x) for x in zip(*sorted(zip(phen_freqs[1:],phen_indices[1:],phen_seqs[1:]),reverse=True))) # exclude undefined phenotype for sorting
phen_ranks = ["x"] + [i+1 for i,x in enumerate(phen_indices[1:])] # rank, largest defined! phenotype has rank 1



# write results ########################################

# write dictionary: genotype index -> NC index to .json file
file_name = results_directory_name + "L.%s_dictionary_genotype_index_to_NC_index.json" % (L)
with open(file_name,'wb') as handle:
    json.dump(gNC_map,handle,sort_keys=True)

# write phenotype characteristics to .txt file
table = [[phen_ranks[i],phen_indices[i],phen_seqs[i],phen_freqs[i]] for i,x in enumerate(phen_ranks)]
file_name = results_directory_name + "L.%s_phenotype_characteristics.txt" % (L)
f = open(file_name,"w")
f.write(tabulate(table,headers=["phen. rank","phen. index","phen. seq.","phen. freq."],floatfmt=(".0f",".0f","",".0f")))
f.close()

# write NC characteristics to .txt file
table = [[NC_ranks[i],NC_indices[i],NC_phen_indices[i],NC_phen_seqs[i],NC_sizes[i],NC_robs[i]] for i,x in enumerate(NC_ranks)]
file_name = results_directory_name + "L.%s_neutral_component_characteristics.txt" % (L)
f = open(file_name,"w")
f.write(tabulate(table,headers=["NC rank","NC index","NC phen. index","NC phen. seq.","NC size","NC rob."],floatfmt=(".0f",".0f",".0f","",".0f")))
f.close()

# write NC average neutral mutations per site to .txt file
table = [[NC_ranks[i],NC_indices[i],NC_phen_indices[i],NC_phen_seqs[i]] for i,x in enumerate(NC_ranks)]
for i,x in enumerate(NC_ranks):
    for l in range(L):
        table[i].append(NC_avg_neutral_mut_per_site_list[i][l])
file_name = results_directory_name + "L.%s_neutral_component_average_neutral_mutations_per_site.txt" % (L)
f = open(file_name,"w")
headers=["NC rank","NC index","NC phen. index","NC phen. seq."] + [""]*L
f.write(tabulate(table,headers=headers))
f.close()

# write NC SD neutral mutations per site to .txt file
table = [[NC_ranks[i],NC_indices[i],NC_phen_indices[i],NC_phen_seqs[i]] for i,x in enumerate(NC_ranks)]
for i,x in enumerate(NC_ranks):
    for l in range(L):
        table[i].append(NC_SD_neutral_mut_per_site_list[i][l])
file_name = results_directory_name + "L.%s_neutral_component_SD_neutral_mutations_per_site.txt" % (L)
f = open(file_name,"w")
headers=["NC rank","NC index","NC phen. index","NC phen. seq."] + [""]*L
f.write(tabulate(table,headers=headers))
f.close()



print "done!"

quit()
