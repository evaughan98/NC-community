import sys
import os
import random
import numpy as np
import json
import networkx as nx



# parameters ########################################

# command line arguments
L = int(sys.argv[1]) # sequence length

# fixed
alphabet = ["A","C","G","U"] # RNA alphabet
random_seed = 1
random.seed(random_seed)
np.random.seed(random_seed)



# directory and file structures ########################################

# input directory name
current_directory_name = os.getcwd()
current_directory_name_oneback = os.path.dirname(current_directory_name)
input_directory_name = current_directory_name_oneback + "/2_NC_analysis/step_3/results/L.%s/" % (L)

# results directory name
results_directory_name = "results"
results_directory_name += "/L.%s" % (L)
results_directory_name += "/"

# create results directory (if not exist)
try:
    os.makedirs(results_directory_name)
except OSError:
    pass



# load input data ########################################

# NC characteristics
input_file_name_1 = input_directory_name + "L.%s_neutral_component_characteristics.txt" % (L)

NC_ranks = list(np.loadtxt(input_file_name_1, usecols=(0,), dtype=int, skiprows=2, unpack=True))
NC_indices = list(np.loadtxt(input_file_name_1, usecols=(1,), dtype=int, skiprows=2, unpack=True))
NC_phen_indices = list(np.loadtxt(input_file_name_1, usecols=(2,), dtype=int, skiprows=2, unpack=True))
NC_phen_seqs = list(np.loadtxt(input_file_name_1, usecols=(3,), dtype=str, skiprows=2, unpack=True))
NC_sizes = list(np.loadtxt(input_file_name_1, usecols=(4,), dtype=int, skiprows=2, unpack=True))

# dictionary: genotype index -> NC index
input_file_name_2 = input_directory_name + "L.%s_dictionary_genotype_index_to_NC_index.json" % (L)
with open(input_file_name_2, 'r') as handle:
    gNC_map = json.load(handle)

# set up NCs = lists of genotype indices
NCs = [[] for NC_index in NC_indices]
for gen_index,NC_index in gNC_map.items():
    i = NC_indices.index(NC_index)
    NCs[i].append(int(gen_index))



# help functions ########################################

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
    
# help function: find neutral neighbour genotype indices for a genotype index
def find_neutral_nb_gen_indices(gen_index,NC,alphabet,L):
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
                if mut_gen_index in NC:
                    neutral_nb_gen_indices.add(mut_gen_index)
    return neutral_nb_gen_indices
    


# program ########################################

# go through all NCs
for index,NC_rank in enumerate(NC_ranks):

    NC = NCs[index]

    # generate graph
    G = nx.Graph()
    for gen_index in NC:
        neutral_nb_gen_indices = find_neutral_nb_gen_indices(gen_index,NC,alphabet,L)
        for neutral_nb_gen_index in neutral_nb_gen_indices:
            G.add_edge(gen_index,neutral_nb_gen_index)
            
    # force-directed graph drawing (via spring layout)
    node_pos = nx.spring_layout(G,iterations=100,scale=2.0)
    
    # add node positions as node attributes
    for node,(x,y) in node_pos.items():
        G.node[node]['pos_x'] = float(x)
        G.node[node]['pos_y'] = float(y)
    
    # save graph
    graph_file_name = results_directory_name + "L.%s_neutral_component_graph_NC.rank.%s.gpickle" % (L,NC_rank)
    nx.write_gpickle(G,graph_file_name)



print "done!"

quit()
