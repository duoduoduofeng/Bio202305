import sys
import calculate_simi
import os

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
        current_block_info_line = None

        for line in file:
            row_id += 1
            line = line.strip()

            # For debug
            if debug_rows_count > 0 and row_id == debug_rows_count:
                break
            
            if not line:
                continue

            # The first block's original information.
            if row_id == 2:
                current_block_info_line = line

            # If starts with "#", new block starts.
            # Thus calculate the average of last block.
            if line.startswith("#"):
                if block_count <= 0:
                    # print(f"Skip line {row_id} which should be the starting.")
                    continue
                
                # parse block info from block info line
                block_info = current_block_info_line.strip().\
                    strip('#').split('\t')
                next_block_info_line = line

                # calculate the average similarity
                # Sometimes the block_info is missing.
                if len(block_info) < 6:
                    print(f"******** WARNING ********\n\
                          file_path: {file_path}, Row_id: {row_id}, current_block_info_line:\n{current_block_info_line}")
                    blocks_avg_sim.append(("NA", "NA", "NA", "NA", "NA", 
                                       current_genome_seq, 
                                       current_genome_seq_2, 
                                       block_count, 
                                       block_sim_sum / block_count))
                else:
                    blocks_avg_sim.append((block_info[0],
                                        block_info[2],
                                        block_info[3],
                                        block_info[4],
                                        block_info[5].split(" ")[0],
                                        current_genome_seq, 
                                        current_genome_seq_2, 
                                        block_count, 
                                        block_sim_sum / block_count))
                    
                # re-init
                current_genome_seq = genome_default
                current_genome_seq_2 = genome_default_2
                block_sim_sum = 0
                block_count = 0
                current_block_info_line = next_block_info_line

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


def save_similarity(output_file, blocks_avg_sim):
    with open(output_file, 'w') as out:
        out.write(f"Block_id\tChromosome_01\tChromosome_02\t\
                  Orientation\tPrior_block_size\t\
                  Genome Sequence\tGenome Sequence 2\t\
                  Block Size\tAverage Percent Identity\n")
        for block in blocks_avg_sim:
            out.write(f"{block[0]}\t{block[1]}\t{block[2]}\t\
                      {block[3]}\t{block[4]}\t\
                        {block[5]}\t{block[6]}\t\
                            {block[7]}\t{block[8]: .2f}\n")
            

def main(args):
    species_list = calculate_simi.parse_species_list(args)
    if species_list is None or len(species_list) <= 0:
        print("No valid input.")
        return

    print(f"The program will deal with the following files: {species_list}")

    # Check the output directory.
    current_dir = os.path.dirname(os.path.abspath(__file__))
    directory = "../triplet_output"
    if not os.path.exists(directory):
        os.makedirs(directory)

    # species_list = ['rainbow_trout', 'chum_salmon']
    for species in species_list:
        input_file = os.path.join(current_dir, species)
        blocks_avg_sim = parse_dagchainer_output(input_file)

        # Print extracted chain information
        output_file = os.path.join(current_dir, directory, 
                                   os.path.basename(species) 
                                   + '.similarity')
        save_similarity(output_file, blocks_avg_sim)


if __name__ == "__main__":
    main(sys.argv)