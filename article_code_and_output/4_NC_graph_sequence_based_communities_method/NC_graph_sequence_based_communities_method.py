import sys
import os
import random
import numpy as np
import networkx as nx
import community
from tabulate import tabulate



# parameters ########################################

# command line arguments
L = int(sys.argv[1]) # sequence length

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
    os.makedirs(results_directory_name + "data/")
except OSError:
    pass
try:
    os.makedirs(results_directory_name + "coarse_grained/graph_files/")
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

# NC average neutral mutations per site
input_file_name_2 = input_directory_name_1 + "L.%s_neutral_component_average_neutral_mutations_per_site.txt" % (L)

NC_average_neutral_mut_per_site_list = [[] for rank in NC_ranks]
for l in range(L):
    load = np.loadtxt(input_file_name_2, usecols=(4+l,), skiprows=2, unpack=True)
    for i,rank in enumerate(NC_ranks):
        NC_average_neutral_mut_per_site_list[i].append(load[i])



# help functions ########################################

# help function (base transfer): genotype index -> genotype sequence
def gen_index_to_gen_seq(gen_index,alphabet,L):
    gen_seq = []
    for site_index in range(L):
        value = gen_index // len(alphabet)**(L-site_index-1)
        letter = alphabet[value]
        gen_seq.append(letter)
        gen_index -= value*len(alphabet)**(L-site_index-1)
    return ''.join(gen_seq)

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

# results lists
results_NC_max_modularity_list = []
results_NC_max_modularity_step_number_list = []

# go through all NCs
for index,NC_rank in enumerate(NC_ranks):


    # exclude NCs with only one node (no edges), for which modularity is not defined
    if NC_sizes[index]==1:
        
        results_NC_max_modularity_list.append("n/a")
        results_NC_max_modularity_step_number_list.append("n/a")


    else:
    
        # load NC graph (and node positions)
        input_file_name_3 = input_directory_name_2 + "L.%s_neutral_component_graph_NC.rank.%s.gpickle" % (L,NC_rank)
        G = nx.read_gpickle(input_file_name_3)
        G_node_pos = {}
        for node in G.nodes():
            G_node_pos[node] = np.array([G.node[node]['pos_x'],G.node[node]['pos_y']])
        
        
        # find constrained positions steps wrt decreasing constraint (increasing average number of neutral mutations)
        site_indices = range(L)
        site_constraints = NC_average_neutral_mut_per_site_list[index]

        sorted_site_constraints, sorted_site_indices = (list(x) for x in zip(*sorted(zip(site_constraints, site_indices))))

        ordered_constrained_positions = [[]] # list that stores for each step the list of corresponding site indices
        check_first_non_fully_constrained = False
        for i,sorted_site_constraint in enumerate(sorted_site_constraints):
            if sorted_site_constraint!=0:
                if check_first_non_fully_constrained==False:
                    ordered_constrained_positions[-1].append(sorted_site_indices[i])
                    check_first_non_fully_constrained = True
                else:
                    # if constraint same as before, add site to same step
                    if sorted_site_constraint==sorted_site_constraints[sorted_site_indices.index(ordered_constrained_positions[-1][-1])]:
                        ordered_constrained_positions[-1].append(sorted_site_indices[i])
                    # if constraint not the same as before, add site to a new step
                    else:
                        ordered_constrained_positions.append([sorted_site_indices[i]])
                        
        # write constrained positions step results to .txt file
        write_step_number = [] # step number = step index + 1, = "x" for fully constrained sites
        for site_index in range(L):
            if site_constraints[site_index]!=0:
                for step_index,constrained_positions in enumerate(ordered_constrained_positions):
                    if site_index in constrained_positions:
                        write_step_number.append(step_index+1)
            else:
                write_step_number.append("x")
        table = [[site_index,NC_average_neutral_mut_per_site_list[index][site_index],write_step_number[site_index]] for site_index in range(L)]
        file_name = results_directory_name + "data/"
        file_name += "L.%s_neutral_component_graph_sequence_based_communities_method_constrained_positions_step_results_NC.rank.%s.txt" % (L,NC_rank)
        f = open(file_name,"w")
        f.write(tabulate(table,headers=["site index","avg. # neutral mut.","step #"]))
        f.close()
        
        
        # go through constrained positions steps, and find for each the respective community structure, modularity and coarse-grained graph
        combined_constrained_positions = []
        step_modularities = []
        check_stop_detection = False
        
        for step_index,constrained_positions in enumerate(ordered_constrained_positions):
        
            # stop if modularity already three times decreased (to reduce computation time)
            if (check_stop_detection==False and (step_index>3 and step_modularities[-1]<step_modularities[-2] and step_modularities[-2]<step_modularities[-3] and step_modularities[-3]<step_modularities[-4])) or check_stop_detection==True:
                
                check_stop_detection = True
                step_modularities.append("n.c.") # not calculated
            
            else:
                # add positions of sites of new step
                combined_constrained_positions.extend(constrained_positions)
                
                # find communities and community sequence building blocks
                community_sequence_building_blocks = [] # list of community sequence building blocks
                community_sequence_building_blocks_seq_counter = [] # list of number of sequences (nodes) belonging to the respective community sequence building block
                part = {}
                
                # go through all nodes
                for node in G.nodes():
                
                    # convert to respective genotype sequence
                    node_sequence = gen_index_to_gen_seq(node,alphabet,L)
                    # replace letters at all "unconstrained" positions by an "x"
                    mod_node_sequence = list(node_sequence)
                    for site_index in range(L):
                        if site_index not in combined_constrained_positions:
                            mod_node_sequence[site_index]="x"
                    mod_node_sequence = ''.join(mod_node_sequence)
                    
                    # check if respective community sequence building block not yet found
                    if mod_node_sequence not in community_sequence_building_blocks:
                        community_sequence_building_blocks.append(mod_node_sequence)
                        community_sequence_building_blocks_seq_counter.append(0)
                        
                    # add node to respective community
                    community_index = community_sequence_building_blocks.index(mod_node_sequence)
                    community_sequence_building_blocks_seq_counter[community_index] += 1
                    part[node] = community_index
                    
                # write community sequence building blocks results to .txt file
                table = [[community_sequence_building_blocks[i],community_sequence_building_blocks_seq_counter[i]] for i,community_sequence_building_block in enumerate(community_sequence_building_blocks)]
                file_name = results_directory_name + "data/"
                file_name += "L.%s_neutral_component_graph_sequence_based_communities_method_community_sequence_building_blocks_NC.rank.%s_step.%s.txt" % (L,NC_rank,step_index+1)
                f = open(file_name,"w")
                f.write(tabulate(table,headers=["community seq. building block","# sequences"]))
                f.close()
                
                # calculate modularity
                modularity = community.modularity(part,G)
                step_modularities.append(modularity)
                
                # generate NC coarse-grained graph
                H = generate_NC_graph_coarse_grained(community_sequence_building_blocks,community_sequence_building_blocks_seq_counter)
                
                # save NC coarse-grained graph file
                graph_file_name = results_directory_name + "coarse_grained/graph_files/" + "L.%s_neutral_component_graph_sequence_based_communities_method_coarse_grained_graph_NC.rank.%s_step.%s.gpickle" % (L,NC_rank,step_index+1)
                nx.write_gpickle(H,graph_file_name)
            
        
    # write modularity step results to .txt file
    table = [[i+1,step_modularities[i]] for i,step_modularity in enumerate(step_modularities)]
    file_name = results_directory_name + "data/"
    file_name += "L.%s_neutral_component_graph_sequence_based_communities_method_modularity_step_results_NC.rank.%s.txt" % (L,NC_rank)
    f = open(file_name,"w")
    f.write(tabulate(table,headers=["# step","modularity"]))
    f.close()
        
    # find step with maximum modularity
    filtered_step_modularities = [x for x in step_modularities if x!="n.c."]
    max_modularity = max(filtered_step_modularities)
    max_modularity_step_number = step_modularities.index(max_modularity) + 1
    
    results_NC_max_modularity_list.append(max_modularity)
    results_NC_max_modularity_step_number_list.append(max_modularity_step_number)
    
    
    
# write results to .txt files ########################################

table = [[NC_ranks[i],NC_indices[i],NC_phen_indices[i],NC_phen_seqs[i],results_NC_max_modularity_list[i],results_NC_max_modularity_step_number_list[i]] for i,x in enumerate(NC_ranks)]
file_name = results_directory_name + "L.%s_neutral_component_graph_sequence_based_communities_method_modularity_results.txt" % (L)
f = open(file_name,"w")
f.write(tabulate(table,headers=["NC rank","NC index","NC phen. index","NC phen. seq.","max. modularity","step #"]))
f.close()



print "done!"

quit()
