import re
import os
import sys
import matplotlib.pyplot as plt
import fit
import fit_gmm
import find_gmm_cutoff
import cutoff_binormal

def parse_dagchainer_output(file_path, debug_rows_count = 0):
    with open(file_path, 'r') as file:
        # global initialization
        row_id = 0
        blocks_avg_sim = []
        genome_default = "CURR"
        genome_default_2 = "CURR_2"

        # block initialization
        current_genome_seq = genome_default
        current_genome_seq_2 = genome_default_2
        block_sim_sum = 0
        block_count = 0
        
        # to record single similarity
        block_id = 0
        each_simi_dict = {}

        for line in file:
            row_id += 1
            line = line.strip()

            # For debug
            if debug_rows_count > 0 and row_id == debug_rows_count:
                break
            
            if not line:
                continue

            # new block single
            if line.startswith("#Ks"):
                block_id += 1

            # If starts with "#", new block starts.
            # Thus calculate the average of last block.
            if line.startswith("#"):
                if block_count <= 0:
                    # print(f"Skip line {row_id} which should be the starting.")
                    continue
                # calculate the average similarity
                blocks_avg_sim.append((current_genome_seq, 
                                       current_genome_seq_2, 
                                       block_count, 
                                       block_sim_sum / block_count))
                    
                # re-init
                current_genome_seq = genome_default
                current_genome_seq_2 = genome_default_2
                block_sim_sum = 0
                block_count = 0

            # check data validation
            fields = line.split('\t')
            if len(fields) != 12:
                if not line.startswith("#"):
                    print(f"Invalid line {row_id} with {len(fields)} fields.")
                continue

            # chr1||start1||stop1||name1||strand1||type1||db_feature_id1||genome_order1||percent_id1
            # chr2||start2||stop2||name2||strand2||type2||db_feature_id2||genome_order2||percent_id2
            chain_info = fields[3].split('||')
            chain_info_2 = fields[7].split('||')
            if len(chain_info) != 9 or len(chain_info_2) != 9:
                print(f"Invalid line {row_id} with fault \
                      columns for the two genome sequence.")
                continue
            # Since the file records paralog, chr1 should be equal to chr2.
            # if chain_info[0] != chain_info_2[0]:
            #     print(f"Invalid line {row_id} with different genome comparison.")
            #     continue
            
            genome_seq = chain_info[0]
            genome_seq_2 = chain_info_2[0]
            # The genome sequence within a block should be the same.
            if current_genome_seq != genome_default \
                and genome_seq != current_genome_seq:
                print(f"Invalid line {row_id} with \
                      different genome sequence within current block")
                continue
            if current_genome_seq_2 != genome_default_2 \
                and genome_seq_2 != current_genome_seq_2:
                print(f"Invalid line {row_id} with \
                      different genome sequence 2 within current block")
                continue

            current_genome_seq = genome_seq
            current_genome_seq_2 = genome_seq_2
            percent_id1 = float(chain_info[-1])
            block_sim_sum += percent_id1
            block_count += 1

            # single
            if block_id not in each_simi_dict:
                each_simi_dict[block_id] = []
            each_simi_dict[block_id].append(percent_id1)
    
    return blocks_avg_sim, each_simi_dict


def save_similarity(output_file, blocks_avg_sim):
    with open(output_file, 'w') as out:
        out.write(f"Genome Sequence\tGenome Sequence 2\t\
                  Block Size\tAverage Percent Identity\n")
        for block in blocks_avg_sim:
            for i in range(block[2]):
                out.write(f"{block[0]}\t{block[1]}\t\
                          {block[2]}\t{block[3]: .2f}\n")

def draw_sim_plot(species, block_sim_data, plot_saving_file):
    data = []
    for block in block_sim_data:
        data.append(float(format(block[3], ".2f")))
    
    plt.clf()
    plt.hist(data, bins=20, edgecolor='black', alpha=0.7)
    plt.xlabel('Values')
    plt.ylabel('Frequency')
    plt.title(f"Similarity Distribution Plot of {species}")
    plt.grid(True)
    plt.savefig(plot_saving_file)
    # plt.show()

def fit_similarity(block_sim_data, fit_file, fit_paras_file):
    data = []
    for block in block_sim_data:
        data.append(float(format(block[3], ".2f")))
    fit.simulate_distribution(data, fit_file, fit_paras_file)

def fit_gmm_similarity(block_sim_data, fit_file, fit_gmm_paras_file):
    data = []
    for block in block_sim_data:
        data.append(float(format(block[3], ".2f")))
    fit_gmm.fit_gmm(data, fit_file, fit_gmm_paras_file)

def find_gmm_cutoffs(block_sim_data, gmm_cutoffs_file):
    data = []
    for block in block_sim_data:
        data.append(float(format(block[3], ".2f")))
    find_gmm_cutoff.find_gmm_cutoff(data, gmm_cutoffs_file)

def find_binormal_cutoffs(block_sim_data, 
                          binormal_cutoffs_parameters, 
                          binormal_cutoffs_plot):
    data = []
    for block in block_sim_data:
        data.append(float(format(block[3], ".2f")))
    return cutoff_binormal.measure_by_mle(data, 
                                   binormal_cutoffs_parameters, 
                                   binormal_cutoffs_plot)

def parse_species_list(args):
    species_list = []
    
    if len(args) < 2:
        print("Please provide a file, "
              "or several files with semicolon as separator, "
              "or folder as an argument.")
    else:
        paths = sys.argv[1:]
        print(f"Your input path: {paths}.\n")
        for path in paths:
            if os.path.isfile(path):
                species_list.append(path)
            elif os.path.isdir(path):
                for root, _, files in os.walk(path):
                    for file in files:
                        file_path = os.path.join(root, file)
                        species_list.append(file_path)
            else:
                filenames = path.split(";")
                for filename in filenames:
                    if os.path.isfile(filename):
                        species_list.append(filename)
                    else:
                        print(f"Invalid file: {filename} "
                              "inside the list you input.")
                        
    return species_list


def judge_side(similarity, cutoff):
    if similarity <= cutoff:
        return "t1"
    else:
        return "t2"


def stat_pair_t1(simi_file, cutoff, t1_stat_file):
    row_id = 0
    t1_count = {}

    with open(simi_file, 'r') as fin, \
        open(t1_stat_file, 'w') as fout:
        for line in fin:
            row_id += 1
            if not line or row_id <= 1:
                continue

            fields = line.strip().split("\t")
            if len(fields) != 4:
                print(f"Invalid line in similarity file, line {row_id}, not 4 fields.")
                continue
            
            similarity = float(fields[3])
            theside = judge_side(similarity, cutoff)
            if theside not in t1_count:
                t1_count[theside] = 0
            t1_count[theside] += 1
        
        sorted_t1_count = {"t1": t1_count["t1"], "t2": t1_count["t2"]}
        fout.write(f"Pair Statistics: {sorted_t1_count}\n")


def stat_pair_between_t1(each_simi_dict, cutoff, another_t1_file):
    t1_count = {}
    # original_count = 0
    # between_count = 0
    with open(another_t1_file, 'w') as fout:
        for block_id in each_simi_dict:
            simis = each_simi_dict[block_id]
            # original_count += len(simis)
            for i in range(1, len(simis) - 1):
                # between_count += 1
                avg_simi = (simis[i-1] + simis[i]) / 2.0
                theside = judge_side(avg_simi, cutoff)
                if theside not in t1_count:
                    t1_count[theside] = 0
                t1_count[theside] += 1
        
        # print(f"original_count = {original_count}, between_count = {between_count}")
        
        sorted_t1_count = {"t1": t1_count["t1"], "t2": t1_count["t2"]}
        fout.write(f"Another Pair Statistics: {sorted_t1_count}\n")


def remove_duplicate(original_file, filtered_file):
    return

def main(args):
    species_list = parse_species_list(args)
    if species_list is None or len(species_list) <= 0:
        print("No valid input.")
        return

    print(f"\n===================  Extract pairs.  ===================")
    print(f"The program will deal with the following files: {species_list}.\n")

    # Check the output directory.
    current_dir = os.path.dirname(os.path.abspath(__file__))
    directory = "../output"
    if not os.path.exists(directory):
        os.makedirs(directory)

    # species_list = ['rainbow_trout', 'chum_salmon']
    for species in species_list:
        print(f"****** Start dealing with {os.path.basename(species)}. ******")

        input_file = os.path.join(current_dir, species)
        blocks_avg_sim, each_simi_dict = parse_dagchainer_output(input_file)

        # Print extracted chain information
        output_file = os.path.join(current_dir, directory, 
                                   os.path.basename(species) 
                                   + '.similarity')
        save_similarity(output_file, blocks_avg_sim)

        # Draw distribution plot
        # plot_file = os.path.join(current_dir, directory, 
        #                          os.path.basename(species) 
        #                          + '.similarity.jpeg')
        # draw_sim_plot(species, blocks_avg_sim, plot_file)

        # Outofuse: Fit the distribution with 110 distributions in scipy lib.
        # fit_file = os.path.join(current_dir, directory, 
        #                          os.path.basename(species) + 
        #                          '.fit.jpeg')
        # fit_paras_file = os.path.join(current_dir, directory, 
        #                          os.path.basename(species) + 
        #                          '.fit.parameters')
        # fit_similarity(blocks_avg_sim, fit_file, fit_paras_file)

        # Fit the distribution by GMM
        # fit_gmm_file = os.path.join(current_dir, directory, 
        #                          os.path.basename(species) + 
        #                          '.fit_gmm.jpeg')
        fit_gmm_file = None
        fit_gmm_paras_file = os.path.join(current_dir, directory, 
                                 os.path.basename(species) + 
                                 '.fit_gmm.parameters')
        fit_gmm_similarity(blocks_avg_sim, fit_gmm_file, fit_gmm_paras_file)

        # Find the cutoff point by GMM
        gmm_cutoff_paras_files = os.path.join(current_dir, directory, 
                                          os.path.basename(species) + 
                                          '.gmm_cutoff.parameters')
        find_gmm_cutoffs(blocks_avg_sim, gmm_cutoff_paras_files)

        # Find the cutoff point for two arbitrary normal distributions
        binormal_cutoff_paras_files = os.path.join(current_dir, directory, 
                                          os.path.basename(species) + 
                                          '.binormal_cutoff.parameters')
        # binormal_cutoff_plot_files = os.path.join(current_dir, directory, 
        #                                   os.path.basename(species) + 
        #                                   '.binormal_cutoff.jpeg')
        binormal_cutoff_plot_files = None
        best_cutoff = find_binormal_cutoffs(blocks_avg_sim, 
                              binormal_cutoff_paras_files, 
                              binormal_cutoff_plot_files)
        
        # Stat the pair t1 & t2
        t1_stat_file = os.path.join(current_dir, directory, 
                                    os.path.basename(species) + 
                                    '.similarity.t1')
        stat_pair_t1(output_file, best_cutoff, t1_stat_file)

        # Calculate another pair t1
        another_t1_file = os.path.join(current_dir, directory, 
                                    os.path.basename(species) + 
                                    '.between.similarity.t1')
        stat_pair_between_t1(each_simi_dict, best_cutoff, another_t1_file)

    print(f"Mission completed. Please check the results in {directory} folder.")

if __name__ == "__main__":
    main(sys.argv)