import sys
import os
import numpy as np
from tabulate import tabulate



# parameters ########################################

# command line arguments
L = int(sys.argv[1]) # sequence length
max_NC_rank = int(sys.argv[2]) # maximum considered NC rank
max_sample_size = int(sys.argv[3]) # maximum considered sample size
number_samples = int(sys.argv[4]) # number of samples considered for each sampling method, sample size, and NC, respectively



# directory and file structures ########################################

# input directory names
current_directory_name = os.getcwd()
current_directory_name_oneback = os.path.dirname(current_directory_name)
input_directory_name_1 = current_directory_name_oneback + "/2_NC_analysis/step_3/results/L.%s/" % (L)
input_directory_name_2 = current_directory_name_oneback + "/4_NC_graph_sequence_based_communities_method/results/L.%s/" % (L)
input_directory_name_3 = current_directory_name_oneback + "/6_NC_sampling/results/L.%s/" % (L)

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
    os.makedirs(results_directory_name + "sampling_average_number_accessed_communities/")
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
input_file_name_2 = input_directory_name_2 + "L.%s_neutral_component_graph_sequence_based_communities_method_modularity_results.txt" % (L)

NC_max_modularity_step_numbers = list(np.loadtxt(input_file_name_2, usecols=(5,), dtype=str, skiprows=2, unpack=True))
NC_max_modularity_step_numbers = [int(x) if x!="n/a" else x for x in NC_max_modularity_step_numbers]

# NC (reference) number of communities
NC_community_numbers = []
for index,NC_rank in enumerate(NC_ranks):
    if NC_sizes[index]>1:
        input_file_name_3 = input_directory_name_2 + "data/"
        input_file_name_3 += "L.%s_neutral_component_graph_sequence_based_communities_method_community_sequence_building_blocks_NC.rank.%s_step.%s.txt" % (L,NC_rank,NC_max_modularity_step_numbers[index])
        reference_community_sequence_building_blocks = np.loadtxt(input_file_name_3, usecols=(0,), dtype=str, skiprows=2, unpack=True)
        try:
            reference_community_sequence_building_blocks = list(reference_community_sequence_building_blocks)
        except TypeError: # only one input line
             reference_community_sequence_building_blocks = [str(reference_community_sequence_building_blocks)]
        NC_community_numbers.append(len(reference_community_sequence_building_blocks))



# program ########################################

# storage lists
NC_RW_average_number_accessed_communities_list_storage = []
NC_RW_SD_number_accessed_communities_list_storage = []

NC_site_scanning_average_number_accessed_communities_list_storage = []
NC_site_scanning_SD_number_accessed_communities_list_storage = []

# go through all NCs (up to and including maximum NC rank and those larger or equal sample size)
for index,NC_rank in enumerate(NC_ranks):

    if NC_ranks[index]<=max_NC_rank and NC_sizes[index]>=max_sample_size:
    
    
        # RW sampling: number accessed communities
        input_file_name_4 = input_directory_name_3 + "sampling_number_accessed_communities/"
        input_file_name_4 += "L.%s_%s_sampling_number_accessed_communities_NC_rank.%s.txt" % (L,"RW",NC_rank)
        results_RW_number_accessed_communities_lists = []
        for n_sample in range(number_samples):
            results_RW_number_accessed_communities_lists.append(list(np.loadtxt(input_file_name_4, usecols=(n_sample+1,), dtype=int, skiprows=2, unpack=True)))
            
        # averaging over samples
        RW_average_number_accessed_communities_list = [0]*max_sample_size
        for i in range(max_sample_size):
            for n_sample in range(number_samples):
                RW_average_number_accessed_communities_list[i] += results_RW_number_accessed_communities_lists[n_sample][i]
        RW_average_number_accessed_communities_list = [float(x)/number_samples for x in RW_average_number_accessed_communities_list]
        
        NC_RW_average_number_accessed_communities_list_storage.append(RW_average_number_accessed_communities_list)
        
        RW_SD_number_accessed_communities_list = [0]*max_sample_size
        for i in range(max_sample_size):
            for n_sample in range(number_samples):
                RW_SD_number_accessed_communities_list[i] += (RW_average_number_accessed_communities_list[i]-results_RW_number_accessed_communities_lists[n_sample][i])**2
        RW_SD_number_accessed_communities_list = [np.sqrt(float(x)/number_samples) for x in RW_SD_number_accessed_communities_list]
        
        NC_RW_SD_number_accessed_communities_list_storage.append(RW_SD_number_accessed_communities_list)
        
        # write results to .txt file
        table = []
        for i in range(max_sample_size):
            table.append([i+1])
            table[i].append(RW_average_number_accessed_communities_list[i])
            table[i].append(RW_SD_number_accessed_communities_list[i])
        file_name = results_directory_name + "sampling_average_number_accessed_communities/"
        file_name += "L.%s_%s_sampling_average_number_accessed_communities_NC_rank.%s.txt" % (L,"RW",NC_rank)
        f = open(file_name,"w")
        headers=["sample size","avg. # acc. comm.","SD # acc. comm."]
        f.write(tabulate(table,headers=headers))
        f.close()
        
        
        # site scanning sampling: number accessed communities
        input_file_name_5 = input_directory_name_3 + "sampling_number_accessed_communities/"
        input_file_name_5 += "L.%s_%s_sampling_number_accessed_communities_NC_rank.%s.txt" % (L,"site_scanning",NC_rank)
        results_site_scanning_number_accessed_communities_lists = []
        for n_sample in range(number_samples):
            results_site_scanning_number_accessed_communities_lists.append(list(np.loadtxt(input_file_name_5, usecols=(n_sample+1,), dtype=int, skiprows=2, unpack=True)))
            
        # averaging over samples
        site_scanning_average_number_accessed_communities_list = [0]*max_sample_size
        for i in range(max_sample_size):
            for n_sample in range(number_samples):
                site_scanning_average_number_accessed_communities_list[i] += results_site_scanning_number_accessed_communities_lists[n_sample][i]
        site_scanning_average_number_accessed_communities_list = [float(x)/number_samples for x in site_scanning_average_number_accessed_communities_list]
        
        NC_site_scanning_average_number_accessed_communities_list_storage.append(site_scanning_average_number_accessed_communities_list)

        site_scanning_SD_number_accessed_communities_list = [0]*max_sample_size
        for i in range(max_sample_size):
            for n_sample in range(number_samples):
                site_scanning_SD_number_accessed_communities_list[i] += (site_scanning_average_number_accessed_communities_list[i]-results_site_scanning_number_accessed_communities_lists[n_sample][i])**2
        site_scanning_SD_number_accessed_communities_list = [np.sqrt(float(x)/number_samples) for x in site_scanning_SD_number_accessed_communities_list]
        
        NC_site_scanning_SD_number_accessed_communities_list_storage.append(site_scanning_SD_number_accessed_communities_list)

        # write results to .txt file
        table = []
        for i in range(max_sample_size):
            table.append([i+1])
            table[i].append(site_scanning_average_number_accessed_communities_list[i])
            table[i].append(site_scanning_SD_number_accessed_communities_list[i])
        file_name = results_directory_name + "sampling_average_number_accessed_communities/"
        file_name += "L.%s_%s_sampling_average_number_accessed_communities_NC_rank.%s.txt" % (L,"site_scanning",NC_rank)
        f = open(file_name,"w")
        headers=["sample size","avg. # acc. comm.","SD # acc. comm."]
        f.write(tabulate(table,headers=headers))
        f.close()

            

# averaging over NCs and write to .txt file

average_NC_RW_average_fraction_accessed_communities_list = [0]*max_sample_size
SD_NC_RW_average_fraction_accessed_communities_list = [0]*max_sample_size

average_NC_site_scanning_average_fraction_accessed_communities_list = [0]*max_sample_size
SD_NC_site_scanning_average_fraction_accessed_communities_list = [0]*max_sample_size

counter = 0
for index,NC_rank in enumerate(NC_ranks):
    if NC_ranks[index]<=max_NC_rank and NC_sizes[index]>=max_sample_size:
        counter += 1
        for i in range(max_sample_size):
            average_NC_RW_average_fraction_accessed_communities_list[i] += float(NC_RW_average_number_accessed_communities_list_storage[index][i])/NC_community_numbers[index]
            average_NC_site_scanning_average_fraction_accessed_communities_list[i] += float(NC_site_scanning_average_number_accessed_communities_list_storage[index][i])/NC_community_numbers[index]
average_NC_RW_average_fraction_accessed_communities_list = [float(x)/counter for x in average_NC_RW_average_fraction_accessed_communities_list]
average_NC_site_scanning_average_fraction_accessed_communities_list = [float(x)/counter for x in average_NC_site_scanning_average_fraction_accessed_communities_list]

counter = 0
for index,NC_rank in enumerate(NC_ranks):
    if NC_ranks[index]<=max_NC_rank and NC_sizes[index]>=max_sample_size:
        counter += 1
        for i in range(max_sample_size):
            SD_NC_RW_average_fraction_accessed_communities_list[i] += (average_NC_RW_average_fraction_accessed_communities_list[i]-float(NC_RW_average_number_accessed_communities_list_storage[index][i])/NC_community_numbers[index])**2
            SD_NC_site_scanning_average_fraction_accessed_communities_list[i] += (average_NC_site_scanning_average_fraction_accessed_communities_list[i]-float(NC_site_scanning_average_number_accessed_communities_list_storage[index][i])/NC_community_numbers[index])**2
SD_NC_RW_average_fraction_accessed_communities_list = [np.sqrt(float(x)/counter) for x in SD_NC_RW_average_fraction_accessed_communities_list]
SD_NC_site_scanning_average_fraction_accessed_communities_list = [np.sqrt(float(x)/counter) for x in SD_NC_site_scanning_average_fraction_accessed_communities_list]


table = []
for i in range(max_sample_size):
    table.append([i+1])
    table[i].append(average_NC_RW_average_fraction_accessed_communities_list[i])
    table[i].append(SD_NC_RW_average_fraction_accessed_communities_list[i])
file_name = results_directory_name
file_name += "L.%s_%s_sampling_average_fraction_accessed_communities_NC_average.txt" % (L,"RW")
f = open(file_name,"w")
headers=["sample size","avg. avg. frac. acc. comm.","SD avg. frac. acc. comm."]
f.write(tabulate(table,headers=headers))
f.close()

table = []
for i in range(max_sample_size):
    table.append([i+1])
    table[i].append(average_NC_site_scanning_average_fraction_accessed_communities_list[i])
    table[i].append(SD_NC_site_scanning_average_fraction_accessed_communities_list[i])
file_name = results_directory_name
file_name += "L.%s_%s_sampling_average_fraction_accessed_communities_NC_average.txt" % (L,"site_scanning")
f = open(file_name,"w")
headers=["sample size","avg. avg. frac. acc. comm.","SD avg. frac. acc. comm."]
f.write(tabulate(table,headers=headers))
f.close()


    
print "done!"

quit()
