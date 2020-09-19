import sys
import os
import random
import numpy as np
import networkx as nx
import community
import igraph as ig
from tabulate import tabulate



# parameters ########################################

# command line arguments
L = int(sys.argv[1]) # sequence length

# fixed
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
    os.makedirs(results_directory_name + "spinglass/graphml_files/")
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

    

# program ########################################

# results lists
results_NC_modularity_Louvain_list = []
results_NC_modularity_spinglass_list = []

# go through all NCs
for index,NC_rank in enumerate(NC_ranks):


    # exclude NCs with only one node (no edges), for which modularity is not defined
    if NC_sizes[index]==1:
        
        results_NC_modularity_Louvain_list.append("n/a")
        results_NC_modularity_spinglass_list.append("n/a")


    else:

        # load NC graph (and node positions)
        input_file_name_2 = input_directory_name_2 + "L.%s_neutral_component_graph_NC.rank.%s.gpickle" % (L,NC_rank)
        G = nx.read_gpickle(input_file_name_2)
        G_node_pos = {}
        for node in G.nodes():
            G_node_pos[node] = np.array([G.node[node]['pos_x'],G.node[node]['pos_y']])
        

        # community detection: Louvain
        NC_part_Louvain = community.best_partition(G) # dictionary that assigns each node (genotype index) a number (starting from 0) corresponding to the community (number) it is a member of
        NC_modularity_Louvain = community.modularity(NC_part_Louvain,G)

        results_NC_modularity_Louvain_list.append(NC_modularity_Louvain)


        # community detection: spinglass

        # convert networkx graph to igraph

        # write networkx graph to .graphml file
        graphml_file_name = results_directory_name + "spinglass/graphml_files/" + "L.%s_neutral_component_graph_NC.rank.%s.graphml" % (L,NC_rank)
        nx.write_graphml(G,graphml_file_name)

        # read in igraph graph
        G_igraph = ig.read(graphml_file_name,format="graphml") # networkx node "names" (i.e. genotype indices) are stored as a "id" attribute to igraph nodes

        # community detection in igraph
        NC_igraph_community_spinglass = G_igraph.community_spinglass()
        NC_igraph_community_membership = NC_igraph_community_spinglass.membership # list that stores for each igraph node a number (starting from 0) corresponding to the community (number) it is a member of

        # transfer found partition back to networkx graph
        NC_part_spinglass = {} # dictionary that assigns each networkx node (genotype index) a number (starting from 0) corresponding to the community (number) it is a member of
        for igraph_node in G_igraph.vs:
            igraph_node_index = igraph_node.index
            netwokx_node_name = int(igraph_node['id'])
            NC_part_spinglass[netwokx_node_name] = NC_igraph_community_membership[igraph_node_index]

        NC_modularity_spinglass = community.modularity(NC_part_spinglass,G)

        results_NC_modularity_spinglass_list.append(NC_modularity_spinglass)



# write results to .txt files ########################################

table = [[NC_ranks[i],NC_indices[i],NC_phen_indices[i],NC_phen_seqs[i],results_NC_modularity_Louvain_list[i],results_NC_modularity_spinglass_list[i]] for i,x in enumerate(NC_ranks)]
file_name = results_directory_name + "L.%s_neutral_component_graph_classical_community_detection_modularity_results.txt" % (L)
f = open(file_name,"w")
f.write(tabulate(table,headers=["NC rank","NC index","NC phen. index","NC phen. seq.","modularity: Louvain","modularity: spinglass"]))
f.close()



print "done!"

quit()
