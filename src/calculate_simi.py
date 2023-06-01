import re
import os
import matplotlib.pyplot as plt

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

        for line in file:
            row_id += 1
            line = line.strip()

            # For debug
            if debug_rows_count > 0 and row_id == debug_rows_count:
                break
            
            if not line:
                continue

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
    
    return blocks_avg_sim

def draw_sim_plot(species, block_sim_data, plot_saving_file):
    data = []
    for block in block_sim_data:
        data.append(float(format(block[3], ".2f")))

    plt.hist(data, bins=20, edgecolor='black', alpha=0.7)
    plt.xlabel('Values')
    plt.ylabel('Frequency')
    plt.title(f"Similarity Distribution Plot of {species}")
    plt.grid(True)
    plt.savefig(plot_saving_file)
    # plt.show()

if __name__ == "__main__":
    species_list = ['rainbow_trout', 'chum_salmon']
    for species in species_list:
        current_dir = os.path.dirname(os.path.abspath(__file__))
        input_file = os.path.join(current_dir, '..', 'data', species)
        blocks_avg_sim = parse_dagchainer_output(input_file)

        # Draw distribution plot
        plot_file = os.path.join(current_dir, '..', 'output', 
                                species + '.similarity.jpeg')
        draw_sim_plot(species, blocks_avg_sim, plot_file)

        # Print extracted chain information
        output_file = os.path.join(current_dir, '..', 'output', 
                                species + '.similarity')
        with open(output_file, 'w') as out:
            out.write(f"Genome Sequence\tGenome Sequence 2\t\
                Block Size\tAverage Percent Identity\n")
            for block in blocks_avg_sim:
                out.write(f"{block[0]}\t{block[1]}\t\
                    {block[2]}\t{block[3]: .2f}\n")