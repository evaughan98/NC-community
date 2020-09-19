import os
import random
import numpy as np
import networkx as nx
from tabulate import tabulate



# parameters ########################################

# fixed
input_fRNAdb_ID = "FR039335"
L = 45

sample_size_list = [1000,2000,3000,4000,5000,6000,7000,8000,9000,10000,20000,30000,40000,50000,60000,70000,80000,90000,100000,200000,300000,400000,500000,600000,700000,800000,900000,1000000]
random_subsample_size_list = [100,1000,"full"]

further_constrained_positions = [14,15,16,22,23,24] # starting from 0

further_constrained_sample_size_list = [10000,100000]
further_constrained_random_subsample_size_list = ["full"]

random_seed = 1
random.seed(random_seed)
np.random.seed(random_seed)

example_sample_size = 100000
example_random_subsample_size = "full"

# alpha (for NC size estimation correction) function parameters
fit_A = 0.68
fit_B = 0.079
alpha = fit_A*(1-np.exp(-fit_B*L))



# directory and file structures ########################################

# input directory name
current_directory_name = os.getcwd()
current_directory_name_oneback = os.path.dirname(current_directory_name)
input_directory_name = current_directory_name_oneback + "/8_fRNAdb_NC_sampling_and_estimations/results/%s/" % (input_fRNAdb_ID)

# results directory name
results_directory_name = "results"
results_directory_name += "/%s" % (input_fRNAdb_ID)
results_directory_name += "/"

# create results directory (if not exist)
try:
    os.makedirs(results_directory_name)
except OSError:
    pass



# help functions ########################################

# help function: generate coarse-grained NC graph
def generate_NC_graph_coarse_grained(community_sequence_building_blocks,community_sequence_building_blocks_seq_counter):

    # generate coarse-grained NC graph
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

# help function: NC size estimation and extrapolated NN size estimation
def NC_size_estimation(ref_structure,sample_avg_neutral_mut_per_site,sample_SD_neutral_mut_per_site,L,alpha):

    NC_size_est = 1.0
    
    for l in range(L):
        if ref_structure[l]==".":
            if (1+sample_avg_neutral_mut_per_site[l]+alpha*sample_SD_neutral_mut_per_site[l])>4:
                NC_size_est = NC_size_est*4.0
            else:
                NC_size_est = NC_size_est*(1+sample_avg_neutral_mut_per_site[l]+alpha*sample_SD_neutral_mut_per_site[l])
        if ref_structure[l]!=".":
            if (1+sample_avg_neutral_mut_per_site[l]+alpha*sample_SD_neutral_mut_per_site[l])>2:
                NC_size_est = NC_size_est*2.0
            else:
                NC_size_est = NC_size_est*(1+sample_avg_neutral_mut_per_site[l]+alpha*sample_SD_neutral_mut_per_site[l])
                
    return NC_size_est
    


# program ########################################

# load input structure
input_file_name_1 = input_directory_name
input_file_name_1 += "%s_fRNAdb_predicted_structure.txt" % (input_fRNAdb_ID)
input_structure = str(np.loadtxt(input_file_name_1, usecols=(0,), dtype=str, skiprows=0, unpack=True))



# 1. # found communities (community sequence building blocks)

# results lists
results_sample_size_list = []
results_random_subsample_size_list = []
results_number_coarse_grained_communities_list = []
results_number_further_coarse_grained_communities_list = []

# go through all sample sizes
for i,sample_size in enumerate(sample_size_list):

    # go through all random subsample sizes
    for j,random_subsample_size in enumerate(random_subsample_size_list):
    
        if random_subsample_size=="full" or sample_size>=random_subsample_size:
         
            # load: sample estimated community sequence building blocks
            input_file_name_1 = input_directory_name
            input_file_name_1 += "%s_fRNAdb_NC_sampling_sample_size.%s_random_subsample_size.%s_sample_estimated_community_sequence_building_blocks.txt" % (input_fRNAdb_ID,sample_size,random_subsample_size)
            community_sequence_building_blocks = list(np.loadtxt(input_file_name_1, usecols=(0,), dtype=str, skiprows=2, unpack=True))
            community_sequence_building_blocks_seq_counter = list(np.loadtxt(input_file_name_1, usecols=(1,), dtype=int, skiprows=2, unpack=True))
            
            number_coarse_grained_communities = len(community_sequence_building_blocks)
             
            # from further constrained positions: find further coarse-grained community sequence building blocks

            # (estimate /) find further coarse-grained communities and community sequence building blocks
            further_coarse_grained_community_sequence_building_blocks = [] # list of further coarse-grained community sequence building blocks
            further_coarse_grained_community_sequence_building_blocks_seq_counter = [] # list of number of sequences (nodes) belonging to the respective further coarse-grained community sequence building block
             
            # go through (initial) community sequence building blocks
            for j,community_sequence_building_block in enumerate(community_sequence_building_blocks):
             
                # replace letters at all "unconstrained" positions by an "x"
                mod_community_sequence_building_block = list(community_sequence_building_block)
                for site_index in range(L):
                    if site_index not in further_constrained_positions:
                        mod_community_sequence_building_block[site_index]="x"
                mod_community_sequence_building_block = ''.join(mod_community_sequence_building_block)
                 
                # check if respective further coarse-grained community sequence building block not yet found
                if mod_community_sequence_building_block not in further_coarse_grained_community_sequence_building_blocks:
                    further_coarse_grained_community_sequence_building_blocks.append(mod_community_sequence_building_block)
                    further_coarse_grained_community_sequence_building_blocks_seq_counter.append(0)
                     
                # add sequences to respective further coarse-grained community
                community_index = further_coarse_grained_community_sequence_building_blocks.index(mod_community_sequence_building_block)
                further_coarse_grained_community_sequence_building_blocks_seq_counter[community_index] += community_sequence_building_blocks_seq_counter[j]
                 
            number_further_coarse_grained_communities = len(further_coarse_grained_community_sequence_building_blocks)
                
            results_sample_size_list.append(sample_size)
            results_random_subsample_size_list.append(random_subsample_size)
            results_number_coarse_grained_communities_list.append(number_coarse_grained_communities)
            results_number_further_coarse_grained_communities_list.append(number_further_coarse_grained_communities)
            
            
# write results to .txt file
table = [[results_sample_size_list[i],results_random_subsample_size_list[i],results_number_coarse_grained_communities_list[i],results_number_further_coarse_grained_communities_list[i]] for i,results_sample_size in enumerate(results_sample_size_list)]
file_name = results_directory_name
file_name += "%s_fRNAdb_NC_sampling_number_sample_estimated_coarse_grained_communities.txt" % (input_fRNAdb_ID)
f = open(file_name,"w")
f.write(tabulate(table,headers=["sample size","random subsample size","# communities","# further coarse-grained communities"]))
f.close()



# 2. further coarse-graining (for specific cases)

# go through specific sample sizes
for i,sample_size in enumerate(further_constrained_sample_size_list):

    # go through specific random subsample sizes
    for j,random_subsample_size in enumerate(further_constrained_random_subsample_size_list):
    
        if random_subsample_size=="full" or sample_size>=random_subsample_size:
        
            # load: sample estimated community sequence building blocks
            input_file_name_2 = input_directory_name
            input_file_name_2 += "%s_fRNAdb_NC_sampling_sample_size.%s_random_subsample_size.%s_sample_estimated_community_sequence_building_blocks.txt" % (input_fRNAdb_ID,sample_size,random_subsample_size)
            community_sequence_building_blocks = list(np.loadtxt(input_file_name_2, usecols=(0,), dtype=str, skiprows=2, unpack=True))
            community_sequence_building_blocks_seq_counter = list(np.loadtxt(input_file_name_2, usecols=(1,), dtype=int, skiprows=2, unpack=True))

            # from further constrained positions: find further coarse-grained community sequence building blocks + graph

            # (estimate /) find further coarse-grained communities and community sequence building blocks
            further_coarse_grained_community_sequence_building_blocks = [] # list of further coarse-grained community sequence building blocks
            further_coarse_grained_community_sequence_building_blocks_seq_counter = [] # list of number of sequences (nodes) belonging to the respective further coarse-grained community sequence building block
            
            # go through (initial) community sequence building blocks
            for j,community_sequence_building_block in enumerate(community_sequence_building_blocks):
            
                # replace letters at all "unconstrained" positions by an "x"
                mod_community_sequence_building_block = list(community_sequence_building_block)
                for site_index in range(L):
                    if site_index not in further_constrained_positions:
                        mod_community_sequence_building_block[site_index]="x"
                mod_community_sequence_building_block = ''.join(mod_community_sequence_building_block)
                
                # check if respective further coarse-grained community sequence building block not yet found
                if mod_community_sequence_building_block not in further_coarse_grained_community_sequence_building_blocks:
                    further_coarse_grained_community_sequence_building_blocks.append(mod_community_sequence_building_block)
                    further_coarse_grained_community_sequence_building_blocks_seq_counter.append(0)
                    
                # add sequences to respective further coarse-grained community
                community_index = further_coarse_grained_community_sequence_building_blocks.index(mod_community_sequence_building_block)
                further_coarse_grained_community_sequence_building_blocks_seq_counter[community_index] += community_sequence_building_blocks_seq_counter[j]
                
                
            # write results to .txt files
            
            # sample estimated further coarse-grained community sequence building blocks
            table = [[further_coarse_grained_community_sequence_building_blocks[k],further_coarse_grained_community_sequence_building_blocks_seq_counter[k]] for k,further_coarse_grained_community_sequence_building_block in enumerate(further_coarse_grained_community_sequence_building_blocks)]
            file_name = results_directory_name
            file_name += "%s_fRNAdb_NC_sampling_sample_size.%s_random_subsample_size.%s_sample_estimated_further_coarse_grained_community_sequence_building_blocks.txt" % (input_fRNAdb_ID,sample_size,random_subsample_size)
            f = open(file_name,"w")
            f.write(tabulate(table,headers=["further coarse-grained community seq. building block","# sequences"]))
            f.close()
        
        
            # generate sample estimated further coarse-grained graph
            H = generate_NC_graph_coarse_grained(further_coarse_grained_community_sequence_building_blocks,further_coarse_grained_community_sequence_building_blocks_seq_counter)

            # save sample estimated further coarse-grained graph file
            graph_file_name = results_directory_name
            graph_file_name += "%s_fRNAdb_NC_sampling_sample_size.%s_random_subsample_size.%s_sample_estimated_further_coarse_grained_coarse_grained_graph.gpickle" % (input_fRNAdb_ID,sample_size,random_subsample_size)
            nx.write_gpickle(H,graph_file_name)
   
   
   
# 3. for example: NC size estimation and fraction of accessed genotypes estimation

# load: random subsample average number of neutral mutations per site, SD number of neutral mutations per site
input_file_name_3 = input_directory_name
input_file_name_3 += "%s_fRNAdb_NC_sampling_sample_size.%s_random_subsample_size.%s_sample_neutral_mutations_per_site.txt" % (input_fRNAdb_ID,example_sample_size,example_random_subsample_size)
random_subsample_avg_neutral_mut_per_site = list(np.loadtxt(input_file_name_3, usecols=(1,), skiprows=2, unpack=True))
random_subsample_SD_neutral_mut_per_site = list(np.loadtxt(input_file_name_3, usecols=(2,), skiprows=2, unpack=True))

# NC size estimation
NC_size_est = NC_size_estimation(input_structure,random_subsample_avg_neutral_mut_per_site,random_subsample_SD_neutral_mut_per_site,L,alpha)

# load: sample estimated community sequence building blocks
input_file_name_4 = input_directory_name
input_file_name_4 += "%s_fRNAdb_NC_sampling_sample_size.%s_random_subsample_size.%s_sample_estimated_community_sequence_building_blocks.txt" % (input_fRNAdb_ID,example_sample_size,example_random_subsample_size)
community_sequence_building_blocks = list(np.loadtxt(input_file_name_4, usecols=(0,), dtype=str, skiprows=2, unpack=True))
community_sequence_building_blocks_seq_counter = list(np.loadtxt(input_file_name_4, usecols=(1,), dtype=int, skiprows=2, unpack=True))

# number of accessed genotypes
number_accessed_genotypes = sum(community_sequence_building_blocks_seq_counter)

# fraction of accessed genotypes
fraction_accessed_genotypes = float(number_accessed_genotypes)/NC_size_est
            
# write results to .txt file
table = [[example_sample_size,example_random_subsample_size,NC_size_est,number_accessed_genotypes,fraction_accessed_genotypes]]
file_name = results_directory_name
file_name += "%s_fRNAdb_NC_sampling_estimated_NC_size_and_fraction_of_accessed_genotypes.txt" % (input_fRNAdb_ID)
f = open(file_name,"w")
f.write(tabulate(table,headers=["sample size","random subsample size","NC size est.","# acc. gen.","frac. acc. gen."]))
f.close()
   
   
        
print "done!"

quit()
