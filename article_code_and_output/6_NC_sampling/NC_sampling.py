import sys
import os
import random
import numpy as np
import networkx as nx
from tabulate import tabulate



# parameters ########################################

# command line arguments
L = int(sys.argv[1]) # sequence length
max_NC_rank = int(sys.argv[2]) # maximum considered NC rank
max_sample_size = int(sys.argv[3]) # maximum considered sample size
number_samples = int(sys.argv[4]) # number of samples considered for each sampling method, sample size, and NC, respectively

# fixed
alphabet = ["A","C","G","U"] # RNA alphabet
random_seed = 1
random.seed(random_seed)
np.random.seed(random_seed)



# directory and file structures ########################################

# input directory names
current_directory_name = os.getcwd()
current_directory_name_oneback = os.path.dirname(current_directory_name)
input_directory_name_1 = current_directory_name_oneback + "/2_NC_analysis/step_3/results/L.%s/" % (L)
input_directory_name_2 = current_directory_name_oneback + "/3_NC_graph_setup_and_layout/results/L.%s/" % (L)
input_directory_name_3 = current_directory_name_oneback + "/4_NC_graph_sequence_based_communities_method/results/L.%s/" % (L)

# results directory name
results_directory_name = "results"
results_directory_name += "/L.%s" % (L)
results_directory_name += "/"

# create results (sub)directories (if not exist)
try:
    os.makedirs(results_directory_name)
except OSError:
    pass
try:
    os.makedirs(results_directory_name + "sampling_number_accessed_communities/")
except OSError:
    pass



# load input data ########################################

# NC characteristics
input_file_name_1 = input_directory_name_1 + "L.%s_neutral_component_characteristics.txt" % (L)

NC_ranks = list(np.loadtxt(input_file_name_1, usecols=(0,), dtype=int, skiprows=2, unpack=True))
NC_indices = list(np.loadtxt(input_file_name_1, usecols=(1,), dtype=int, skiprows=2, unpack=True))
NC_phen_indices = list(np.loadtxt(input_file_name_1, usecols=(2,), dtype=int, skiprows=2, unpack=True))
NC_phen_seqs = list(np.loadtxt(input_file_name_1, usecols=(3,), dtype=str, skiprows=2, unpack=True))
NC_sizes = list(np.loadtxt(input_file_name_1, usecols=(4,), dtype=int, skiprows=2, unpack=True))

# NC sequence based communities method results
input_file_name_2 = input_directory_name_3 + "L.%s_neutral_component_graph_sequence_based_communities_method_modularity_results.txt" % (L)

NC_max_modularity_step_numbers = list(np.loadtxt(input_file_name_2, usecols=(5,), dtype=str, skiprows=2, unpack=True))
NC_max_modularity_step_numbers = [int(x) if x!="n/a" else x for x in NC_max_modularity_step_numbers]



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

# help function: 'RW' neutral mutation for a genotype index
def RW_neutral_mutation(gen_index,ref_NC,alphabet,L):
    gen_seq = gen_index_to_gen_seq(gen_index,alphabet,L)
    gen_seq_array = list(gen_seq)
    success = False
    while success==False:
        # random position
        mut_position = random.randint(0,L-1)
        # mutation alphabet
        mut_alphabet = [letter for letter in alphabet if letter!=gen_seq_array[mut_position]]
        # random letter
        mut_letter = random.choice(mut_alphabet)
        mut_gen_seq_array = [j for j in gen_seq_array] # to preserve "gen_seq_array"
        # mutate
        mut_gen_seq_array[mut_position] = mut_letter
        mut_gen_index = gen_seq_to_gen_index(''.join(mut_gen_seq_array),alphabet,L)
        # check if neutral
        if mut_gen_index in ref_NC:
            success=True
    return mut_gen_index

# help function: 'site scanning' neutral mutation for a genotype index
def site_scanning_neutral_mutation(gen_index,ref_mut_position,ref_NC,alphabet,L):
    gen_seq = gen_index_to_gen_seq(gen_index,alphabet,L)
    gen_seq_array = list(gen_seq)
    success = False
    while success==False:
        # position
        mut_position = ref_mut_position
        # mutation alphabet
        mut_alphabet = [letter for letter in alphabet if letter!=gen_seq_array[mut_position]]
        # randomly test mutation alphabet until success
        while len(mut_alphabet)>0:
            # random letter
            mut_letter = random.choice(mut_alphabet)
            mut_gen_seq_array = [j for j in gen_seq_array] # to preserve "gen_seq_array"
            # mutate
            mut_gen_seq_array[mut_position] = mut_letter
            mut_gen_index = gen_seq_to_gen_index(''.join(mut_gen_seq_array),alphabet,L)
            # check if neutral
            if mut_gen_index in ref_NC:
                success=True
                break
            # if no success, update mutation alphabet
            mut_alphabet.remove(mut_letter)
        # if no success, go to next position
        if len(mut_alphabet)==0 and success==False:
            ref_mut_position = (ref_mut_position+1) % L
    return mut_gen_index, mut_position



# program ########################################

# go through all NCs (up to and including maximum NC rank and those larger or equal sample size)
for index,NC_rank in enumerate(NC_ranks):

    if NC_ranks[index]<=max_NC_rank and NC_sizes[index]>=max_sample_size:
    
    
        # load NC graph
        input_file_name_3 = input_directory_name_2 + "L.%s_neutral_component_graph_NC.rank.%s.gpickle" % (L,NC_rank)
        G = nx.read_gpickle(input_file_name_3)
        
        NC = list(G.nodes())
        
        # load reference community sequence building blocks
        input_file_name_4 = input_directory_name_3 + "data/"
        input_file_name_4 += "L.%s_neutral_component_graph_sequence_based_communities_method_community_sequence_building_blocks_NC.rank.%s_step.%s.txt" % (L,NC_rank,NC_max_modularity_step_numbers[index])
        reference_community_sequence_building_blocks = np.loadtxt(input_file_name_4, usecols=(0,), dtype=str, skiprows=2, unpack=True)
        try:
            reference_community_sequence_building_blocks = list(reference_community_sequence_building_blocks)
        except TypeError: # only one input line
             reference_community_sequence_building_blocks = [str(reference_community_sequence_building_blocks)]

        
        # results lists
        results_RW_number_accessed_communities_lists = []
        results_site_scanning_number_accessed_communities_lists = []
        
        
        # go through all samples
        for n_sample in range(number_samples):
        
        
            # 1. sampling and # accessed communities
            
            # RW sampling
            ref_gen_index = random.sample(NC,1)[0]
            RW_sample = [ref_gen_index]
            while len(RW_sample) < max_sample_size:
                mut_gen_index = RW_neutral_mutation(ref_gen_index,NC,alphabet,L)
                RW_sample.append(mut_gen_index)
                ref_gen_index = mut_gen_index
                
            RW_number_accessed_communities = []
            RW_accessed_communities = set([])
            for gen_index in RW_sample:
                gen_sequence = gen_index_to_gen_seq(gen_index,alphabet,L)
                for reference_community_sequence_building_block in reference_community_sequence_building_blocks:
                    check=True
                    for site_index in range(L):
                        if reference_community_sequence_building_block[site_index]!="x" and reference_community_sequence_building_block[site_index]!=gen_sequence[site_index]:
                            check=False
                    if check==True:
                        RW_accessed_communities.add(reference_community_sequence_building_block)
                RW_number_accessed_communities.append(len(RW_accessed_communities))
                
            results_RW_number_accessed_communities_lists.append(RW_number_accessed_communities)
            
            
            # site scanning sampling
            ref_gen_index = random.sample(NC,1)[0]
            site_scanning_sample = [ref_gen_index]
            ref_mut_position = 0
            while len(site_scanning_sample) < max_sample_size:
                mut_gen_index, mut_position = site_scanning_neutral_mutation(ref_gen_index,ref_mut_position,NC,alphabet,L)
                site_scanning_sample.append(mut_gen_index)
                ref_gen_index = mut_gen_index
                ref_mut_position = (mut_position+1) % L
                
            site_scanning_number_accessed_communities = []
            site_scanning_accessed_communities = set([])
            for gen_index in site_scanning_sample:
                gen_sequence = gen_index_to_gen_seq(gen_index,alphabet,L)
                for reference_community_sequence_building_block in reference_community_sequence_building_blocks:
                    check=True
                    for site_index in range(L):
                        if reference_community_sequence_building_block[site_index]!="x" and reference_community_sequence_building_block[site_index]!=gen_sequence[site_index]:
                            check=False
                    if check==True:
                        site_scanning_accessed_communities.add(reference_community_sequence_building_block)
                site_scanning_number_accessed_communities.append(len(site_scanning_accessed_communities))
                
            results_site_scanning_number_accessed_communities_lists.append(site_scanning_number_accessed_communities)
                        
        
            # write example results to .txt file (for first two samples per NC)
            if n_sample==0 or n_sample==1:
                        
                try:
                    os.makedirs(results_directory_name + "examples/NC_rank.%s/sample.%s/" % (NC_rank,n_sample+1))
                except OSError:
                    pass
                            
                            
                # RW sample genotypes
                table = [[k,RW_sample[k],gen_index_to_gen_seq(RW_sample[k],alphabet,L)] for k in range(max_sample_size)]
                file_name = results_directory_name + "examples/NC_rank.%s/sample.%s/" % (NC_rank,n_sample+1)
                file_name += "L.%s_%s_sampling_sample_genotypes_NC_rank.%s_sample.%s.txt" % (L,"RW",NC_rank,n_sample+1)
                f = open(file_name,"w")
                headers=["#","genotype index","genotype sequence"]
                f.write(tabulate(table,headers=headers))
                f.close()
            
                # site scanning sample genotypes
                table = [[k,site_scanning_sample[k],gen_index_to_gen_seq(site_scanning_sample[k],alphabet,L)] for k in range(max_sample_size)]
                file_name = results_directory_name + "examples/NC_rank.%s/sample.%s/" % (NC_rank,n_sample+1)
                file_name += "L.%s_%s_sampling_sample_genotypes_NC_rank.%s_sample.%s.txt" % (L,"site_scanning",NC_rank,n_sample+1)
                f = open(file_name,"w")
                headers=["#","genotype index","genotype sequence"]
                f.write(tabulate(table,headers=headers))
                f.close()
                        
                        
                        
        # write results to .txt file
        
        # RW sampling: number accessed communities
        table = []
        for i in range(max_sample_size):
            table.append([i+1])
            for n_sample in range(number_samples):
                table[i].append(results_RW_number_accessed_communities_lists[n_sample][i])
        file_name = results_directory_name + "sampling_number_accessed_communities/"
        file_name += "L.%s_%s_sampling_number_accessed_communities_NC_rank.%s.txt" % (L,"RW",NC_rank)
        f = open(file_name,"w")
        headers=["sample size"] + ["#%s" % (n_sample+1) for n_sample in range(number_samples)]
        f.write(tabulate(table,headers=headers))
        f.close()

        # site scanning sampling: number accessed communities
        table = []
        for i in range(max_sample_size):
            table.append([i+1])
            for n_sample in range(number_samples):
                table[i].append(results_site_scanning_number_accessed_communities_lists[n_sample][i])
        file_name = results_directory_name + "sampling_number_accessed_communities/"
        file_name += "L.%s_%s_sampling_number_accessed_communities_NC_rank.%s.txt" % (L,"site_scanning",NC_rank)
        f = open(file_name,"w")
        headers=["sample size"] + ["#%s" % (n_sample+1) for n_sample in range(number_samples)]
        f.write(tabulate(table,headers=headers))
        f.close()



print "done!"

quit()
