import sys
import os
import random
import numpy as np
import networkx as nx
from tabulate import tabulate

import RNA



# parameters ########################################

# command line arguments
input_fRNAdb_ID = str(sys.argv[1]) # input fRNAdb ID
input_sequence = str(sys.argv[2]) # input sequence
sample_size_list = [int(sample_size) for sample_size in sys.argv[3][1:-1].split(",")] # list of considered sample sizes
random_subsample_size_list = [int(sample_size) if sample_size!="full" else sample_size for sample_size in sys.argv[4][1:-1].split(",")] # list of considered random subsample sizes (not for random sampling)

max_sample_size = max(sample_size_list)

# fixed
alphabet = ["A","C","G","U"] # RNA alphabet
alphabet_compatible_bp_partners = [["U"],["G"],["C","U"],["A","G"]] # list of (lists of) letters that in principle are able (/ compatible) to form a bp with letter at respective position in RNA alphabet
random_seed = 1
random.seed(random_seed)
np.random.seed(random_seed)



# directory and file structures ########################################

# results directory name
results_directory_name = "results"
results_directory_name += "/%s" % (input_fRNAdb_ID)
results_directory_name += "/"

# create results directory (if not exist)
try:
    os.makedirs(results_directory_name)
except OSError:
    pass



# process + check input ########################################

# sequence length
L = len(input_sequence)

# process input sequence: replace T in sequence by U
sequence_processed = ""
for character in input_sequence:
    if character=="T":
        sequence_processed += "U"
    else:
        sequence_processed += character
input_sequence = sequence_processed

# check if input sequence leads to defined phenotype
(structure, mfe) = RNA.fold(input_sequence)
if structure=="."*L:
    print "error: input sequence leads to undefined phenotype (unbound structure)"
    quit()
ref_structure = structure



# help functions ########################################

# help function: find structure
def find_structure(sequence):
    (structure, mfe) = RNA.fold(sequence)
    return structure

# help function: find bp partner positions
def find_bp_partner_positions(structure):
    bp_partner_positions = ["x" for site in structure] # list that has at each position the position index of the bp partner, "x" if position is not paired
    for i,symbol_1 in enumerate(structure):
        if symbol_1=='(':
            open_counter = 1
            close_counter = 0
            for j,symbol_2 in enumerate(structure):
                if j>i and symbol_2=='(':
                    open_counter += 1
                if j>i and symbol_2==')':
                    close_counter += 1
                if close_counter==open_counter:
                    bp_partner_positions[i]=j
                    bp_partner_positions[j]=i
                    break
    return bp_partner_positions

# help function: 'accelerated site scanning' neutral mutation for a sequence
def site_scanning_neutral_mutation(sequence,ref_mut_position,ref_structure,bp_partner_positions,L,alphabet,alphabet_compatible_bp_partners):
    success = False
    while success==False:
        # position
        mut_position = ref_mut_position
        # mutation alphabet
        mut_alphabet = [letter for letter in alphabet if letter!=sequence[mut_position]]
        # randomly test mutation alphabet until success
        while len(mut_alphabet)>0:
            # random letter
            mut_letter = random.choice(mut_alphabet)
            # check if mutation affects unpaired site or if it affects paired site if it would lead to compatible bp
            if ref_structure[mut_position]=="." or (ref_structure[mut_position]!="." and mut_letter in alphabet_compatible_bp_partners[alphabet.index(sequence[bp_partner_positions[mut_position]])]):
                mut_sequence = [j for j in sequence] # to preserve "sequence"
                # mutate
                mut_sequence[mut_position] = mut_letter
                mut_sequence = "".join(mut_sequence)
                # check if neutral
                mut_structure = find_structure(mut_sequence)
                if mut_structure==ref_structure:
                    success=True
                    break
            # if no success, update mutation alphabet
            mut_alphabet.remove(mut_letter)
        # if no success, go to next position
        if len(mut_alphabet)==0 and success==False:
            ref_mut_position = (ref_mut_position+1) % L
    return mut_sequence, mut_position

# help function: find neutral mutations per site, and neutral neighbour sequences for a sequence
def find_neutral_mut_per_site_and_neutral_nb_sequences(sequence,ref_structure,bp_partner_positions,L,alphabet,alphabet_compatible_bp_partners):
    seq_neutral_mut_per_site = [0]*L
    seq_neutral_nb_sequences = []
    # go through all characters of the sequence
    for i,character in enumerate(sequence):
        # go through all letters of the alphabet
        for letter in alphabet:
            # mutate character to "value" of letter if character != letter
            if character!=letter:
                # check if mutation affects unpaired site or if it affects paired site if it would lead to compatible bp
                if ref_structure[i]=="." or (ref_structure[i]!="." and letter in alphabet_compatible_bp_partners[alphabet.index(sequence[bp_partner_positions[i]])]):
                    mut_sequence = [j for j in sequence] # to preserve "sequence"
                    # mutate
                    mut_sequence[i] = letter
                    mut_sequence = "".join(mut_sequence)
                    # check if neutral
                    mut_structure = find_structure(mut_sequence)
                    if mut_structure==ref_structure:
                        # if neutral mutation, add it to the specific site
                        seq_neutral_mut_per_site[i] += 1
                        seq_neutral_nb_sequences.append(mut_sequence)
    return seq_neutral_mut_per_site, seq_neutral_nb_sequences

# help function: sample measurement (averaging)
def sample_measurement(sample_seq_neutral_mut_per_site_storage,L):

    sample_avg_neutral_mut_per_site = [0 for l in range(L)]
    sample_SD_neutral_mut_per_site = [0 for l in range(L)]

    for k in range(len(sample_seq_neutral_mut_per_site_storage)):
        for l in range(L):
            sample_avg_neutral_mut_per_site[l] += sample_seq_neutral_mut_per_site_storage[k][l]
    for l in range(L):
        sample_avg_neutral_mut_per_site[l] = float(sample_avg_neutral_mut_per_site[l])/len(sample_seq_neutral_mut_per_site_storage)

    for k in range(len(sample_seq_neutral_mut_per_site_storage)):
        for l in range(L):
            sample_SD_neutral_mut_per_site[l] += (sample_seq_neutral_mut_per_site_storage[k][l]-sample_avg_neutral_mut_per_site[l])**2
    for l in range(L):
        if len(sample_seq_neutral_mut_per_site_storage)>1:
            sample_SD_neutral_mut_per_site[l] = np.sqrt(float(sample_SD_neutral_mut_per_site[l])/(len(sample_seq_neutral_mut_per_site_storage)-1))
        if len(sample_seq_neutral_mut_per_site_storage)==1:
            sample_SD_neutral_mut_per_site[l] = 0

    return sample_avg_neutral_mut_per_site, sample_SD_neutral_mut_per_site

# help function: generate NC coarse-grained graph
def generate_NC_graph_coarse_grained(community_sequence_building_blocks,community_sequence_building_blocks_seq_counter):

    # generate NC coarse-grained graph
    H = nx.Graph()

    for index_1,community_sequence_building_block_1 in enumerate(community_sequence_building_blocks):
        H.add_node(index_1)

    for index_1,community_sequence_building_block_1 in enumerate(community_sequence_building_blocks):
        for index_2,community_sequence_building_block_2 in enumerate(community_sequence_building_blocks):
            if index_1!=index_2:
                # if community sequence building blocks only differ by one letter (i.e. by 1-point mutation) connect nodes
                differ_counter = 0
                for site_index in range(L):
                    if community_sequence_building_block_1[site_index]!=community_sequence_building_block_2[site_index]:
                        differ_counter += 1
                if differ_counter==1:
                    H.add_edge(index_1, index_2)

    # relative node sizes
    for node in H.nodes():
        relative_size = float(community_sequence_building_blocks_seq_counter[node])/sum(community_sequence_building_blocks_seq_counter)
        # add as node attribute
        H.node[node]['rel_size'] = relative_size
    
    # node labels
    for node in H.nodes():
        node_label = community_sequence_building_blocks[node]
        # add as node attribute
        H.node[node]['label'] = node_label
    
    return H



# program ########################################

# write predicted structure to .txt file
file_name = results_directory_name
file_name += "%s_fRNAdb_predicted_structure.txt" % (input_fRNAdb_ID)
f = open(file_name,"w")
f.write("%s" % (ref_structure))
f.close()

bp_partner_positions = find_bp_partner_positions(ref_structure)

# sampling (accelerated site scanning) from input sequence
ref_sequence = input_sequence
sample = [ref_sequence]
ref_mut_position = 0
while len(sample) < max_sample_size:
    mut_sequence, mut_position = site_scanning_neutral_mutation(ref_sequence,ref_mut_position,ref_structure,bp_partner_positions,L,alphabet,alphabet_compatible_bp_partners)
    sample.append(mut_sequence)
    ref_sequence = mut_sequence
    ref_mut_position = (mut_position+1) % L

# sequence measurements
sample_seq_neutral_mut_per_site_storage = []
sample_seq_neutral_nb_sequences_storage = []
for sequence in sample:
    seq_neutral_mut_per_site, seq_neutral_nb_sequences = find_neutral_mut_per_site_and_neutral_nb_sequences(sequence,ref_structure,bp_partner_positions,L,alphabet,alphabet_compatible_bp_partners)
    sample_seq_neutral_mut_per_site_storage.append(seq_neutral_mut_per_site)
    sample_seq_neutral_nb_sequences_storage.append(seq_neutral_nb_sequences)
    
# go through all sample sizes
for i,sample_size in enumerate(sample_size_list):

    # go through all random subsample sizes
    for j,random_subsample_size in enumerate(random_subsample_size_list):
    
        if random_subsample_size=="full" or sample_size>=random_subsample_size:
        
            # random subsampling
            if random_subsample_size=="full":
                random_subsample_indices = range(sample_size)
            else:
                random_subsample_indices = random.sample(range(sample_size),random_subsample_size)
            
            random_subsample_seq_neutral_mut_per_site_storage = []
            for k in random_subsample_indices:
                random_subsample_seq_neutral_mut_per_site_storage.append(sample_seq_neutral_mut_per_site_storage[k])
        
            # random subsample measurement
            random_subsample_avg_neutral_mut_per_site, random_subsample_SD_neutral_mut_per_site = sample_measurement(random_subsample_seq_neutral_mut_per_site_storage,L)
            
            # find constrained positions
            constrained_positions = []
            
            for site_index,symbol in enumerate(ref_structure):
                if symbol!=".":
                    constrained_positions.append(site_index)
            
            # determine available sequence pool
            sequence_pool = set([])
            for k in range(sample_size):
                sequence_pool.add(sample[k])
            for k in random_subsample_indices:
                for neutral_nb_sequence in sample_seq_neutral_nb_sequences_storage[k]:
                    sequence_pool.add(neutral_nb_sequence)
            
            
            # (estimate /) find communities and community sequence building blocks
            community_sequence_building_blocks = [] # list of community sequence building blocks
            community_sequence_building_blocks_seq_counter = [] # list of number of sequences (nodes) belonging to the respective community sequence building block
            
            # go through sequences
            for sequence in sequence_pool:
            
                # replace letters at all "unconstrained" positions by an "x"
                mod_sequence = list(sequence)
                for site_index in range(L):
                    if site_index not in constrained_positions:
                        mod_sequence[site_index]="x"
                mod_sequence = ''.join(mod_sequence)
                
                # check if respective community sequence building block not yet found
                if mod_sequence not in community_sequence_building_blocks:
                    community_sequence_building_blocks.append(mod_sequence)
                    community_sequence_building_blocks_seq_counter.append(0)
                    
                # add sequences to respective community
                community_index = community_sequence_building_blocks.index(mod_sequence)
                community_sequence_building_blocks_seq_counter[community_index] += 1
                
                
            # write results to .txt files
            
            # sample neutral mutations per site
            table = []
            for site_index in range(L):
                if site_index in constrained_positions:
                    table.append([site_index,random_subsample_avg_neutral_mut_per_site[site_index],random_subsample_SD_neutral_mut_per_site[site_index],"y"])
                else:
                    table.append([site_index,random_subsample_avg_neutral_mut_per_site[site_index],random_subsample_SD_neutral_mut_per_site[site_index],"n"])
            file_name = results_directory_name
            file_name += "%s_fRNAdb_NC_sampling_sample_size.%s_random_subsample_size.%s_sample_neutral_mutations_per_site.txt" % (input_fRNAdb_ID,sample_size,random_subsample_size)
            f = open(file_name,"w")
            f.write(tabulate(table,headers=["site index","avg. # neutral mut.","SD # neutral mut.","constrained?"]))
            f.close()
            
            # sample estimated community sequence building blocks
            table = [[community_sequence_building_blocks[k],community_sequence_building_blocks_seq_counter[k]] for k,community_sequence_building_block in enumerate(community_sequence_building_blocks)]
            file_name = results_directory_name
            file_name += "%s_fRNAdb_NC_sampling_sample_size.%s_random_subsample_size.%s_sample_estimated_community_sequence_building_blocks.txt" % (input_fRNAdb_ID,sample_size,random_subsample_size)
            f = open(file_name,"w")
            f.write(tabulate(table,headers=["community seq. building block","# sequences"]))
            f.close()
        
        
            # generate sample estimated coarse-grained graph
            H = generate_NC_graph_coarse_grained(community_sequence_building_blocks,community_sequence_building_blocks_seq_counter)
        
            # save sample estimated coarse-grained graph file
            graph_file_name = results_directory_name
            graph_file_name += "%s_fRNAdb_NC_sampling_sample_size.%s_random_subsample_size.%s_sample_estimated_coarse_grained_graph.gpickle" % (input_fRNAdb_ID,sample_size,random_subsample_size)
            nx.write_gpickle(H,graph_file_name)
        


print "done!"

quit()
