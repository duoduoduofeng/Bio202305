import os
import sys
import calculate_simi


def calculate_common_ratio(list1, list2):
    # Convert both lists to sets and find their intersection
    common_elements = set(list1) & set(list2)
    
    # Count the number of common elements
    common_count = len(common_elements)
    
    # Calculate the ratio of common elements to the length of list1
    ratio = common_count / len(list1) if len(list1) > 0 else 0
    return ratio


def remove_dag_dup(dag_input_file, dag_output_file, 
                   dup_threshold = 0.9, debug_rows_count = 0):
    
    with open(dag_input_file, 'r') as fin, \
        open(dag_output_file, 'w') as fout:
        
        # global initialization
        row_id = 0
        genome_default = "CURR"
        genome_default_2 = "CURR_2"

        # block initialization
        current_genome_seq = genome_default
        current_genome_seq_2 = genome_default_2
        block_seqs_dict = {}
        block_count = 0
        block_id = 1
        block_fids = []
        
        # store line batch
        curlines = []

        for line in fin:
            row_id += 1
            line = line.strip()

            # For debug
            if debug_rows_count > 0 and row_id == debug_rows_count:
                break
            
            if not line:
                continue
            curlines.append(line)

            # If starts with "#", new block starts.
            if line.startswith("#Ks"):
                if block_count <= 0:
                    print(f"Skip line {row_id} which should be the starting.")
                    continue
                
                genomes = (current_genome_seq, current_genome_seq_2)
                if genomes not in block_seqs_dict:
                    block_seqs_dict[genomes] = {}
                if block_id not in block_seqs_dict[genomes]:
                    block_seqs_dict[genomes][block_id] = block_fids
                    for aline in curlines:
                        fout.write(f"{aline}\n")
                else:
                    # Judge the duplication.
                    old_fids = block_seqs_dict[genomes][block_id]
                    dup_ratio = calculate_common_ratio(block_fids, old_fids)
                    # print(f"dup_ratio = {dup_ratio}")
                    if dup_ratio < dup_threshold:
                        for aline in curlines:
                            fout.write(f"{aline}\n")
                    else:
                        print(f"Duplicates in {genomes}, block_id = {block_id}, ratio = {dup_ratio}")
    
                # re-init
                current_genome_seq = genome_default
                current_genome_seq_2 = genome_default_2
                block_fids = []
                block_count = 0
                curlines = []

                # For the next block
                block_id += 1
            
            # new genome pairs signal
            # The block id will shift to 1 for the last block in a genome pair,
            # which doesn't affect the result, thus ignore it.
            if line.startswith("#1"):
                block_id = 1

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
            
            genome_seq = chain_info[0]
            genome_seq_2 = chain_info_2[0]
            if genome_seq > genome_seq_2: # rank by alphabet order
                tmp = genome_seq
                genome_seq = genome_seq_2
                genome_seq_2 = tmp
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

            fid_1 = chain_info[6]
            fid_2 = chain_info_2[6]
            if fid_1 <= fid_2:
                block_fids.append((fid_1, fid_2))
            else:
                block_fids.append((fid_2, fid_1))

            block_count += 1
        
        if len(curlines) > 0:
            for aline in curlines:
                fout.write(f"{aline}\n")


def main(args):
    species_list = calculate_simi.parse_species_list(args)
    if species_list is None or len(species_list) <= 0:
        print("No valid input.")
        return

    print(f"\n===================  Extract pairs.  ===================")
    print(f"The program will deal with the following files: {species_list}.\n")

    # Check the output directory.
    current_dir = os.path.dirname(os.path.abspath(__file__))
    directory = "../data/paralogs_outputs"
    if not os.path.exists(directory):
        os.makedirs(directory)

    # species_list = ['rainbow_trout', 'chum_salmon']
    for species in species_list:
        print(f"****** Start dealing with {os.path.basename(species)}. ******")

        input_file = os.path.join(current_dir, species)
        output_file = os.path.join(current_dir, directory, 
                                   os.path.basename(species))
        remove_dag_dup(input_file, output_file)


if __name__ == "__main__":
    main(sys.argv)