import os
import sys
import calculate_simi
import cutoff_binormal


def find_binormal_cutoffs(block_sim_data, 
                          binormal_cutoffs_parameters, 
                          binormal_cutoffs_plot):
    data = []
    for block in block_sim_data:
        data.append(float(format(block[3], ".2f")))
    cutoff = cutoff_binormal.measure_by_mle(data, 
                                   binormal_cutoffs_parameters, 
                                   binormal_cutoffs_plot)
    return cutoff
    

def traverse_each_species(file_path, debug_rows_count = 0):
    genome_pair_sim = {}
    with open(file_path, 'r') as file:
        row_id = 0

        for line in file:
            row_id += 1
            line = line.strip()

            # For debug
            if debug_rows_count > 0 and row_id == debug_rows_count:
                break

            if line.startswith("#"):
                continue

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

            genome_segment = int(chain_info[6])
            genome_segment_2 = int(chain_info_2[6])
            # ensure genome_segment is not greater than genome_segment_2
            if genome_segment > genome_segment_2:
                tmp = genome_segment
                genome_segment = genome_segment_2
                genome_segment_2 = tmp
            similarity = float(chain_info[-1])

            new_pair = (genome_segment, genome_segment_2)
            if new_pair in genome_pair_sim:
                print(f"Invalid line {row_id} which contains the same \
                      genome segment pairs with previous lines.")
                continue
            genome_pair_sim[new_pair] = similarity

    return genome_pair_sim


def judge_side(cutoff, similarity):
    if similarity <= cutoff:
        return "t1"
    else:
        return "t2"


def extract_triples(genome_pair_sim, similarity_cutoff, triple_output_file):
    triples = {}
    pairs = genome_pair_sim.keys()
    # print(pairs)

    #Identify unique elements in the list of pairs
    unique_elements = set()
    for pair in pairs:
        unique_elements.update(pair)

    for pair1 in pairs:
        for pair2 in pairs:
            if pair1 != pair2:
                common_element = None
                for element in pair1:
                    if element in pair2 and element in unique_elements:
                        common_element = element
                        break

                if common_element:
                    remaining_elements = [e for e in pair1 + pair2 if e != common_element]
                    # (A, B), (A, C) and (B, C) should all appear in the pairs
                    if (remaining_elements[0], remaining_elements[1]) in pairs:
                        triple = [common_element] + remaining_elements
                        triple.sort()  # Sort for consistency
                        triple = tuple(triple)
                        if triple not in triples:
                            triples[triple] = [genome_pair_sim[(triple[0], triple[1])], 
                                               genome_pair_sim[(triple[0], triple[2])],
                                               genome_pair_sim[(triple[1], triple[2])]]

    if triple_output_file:
        with open(triple_output_file, 'w') as pf:
            pf.write(f"Cutoff point: {similarity_cutoff}\n")
            for triple in triples:
                sims = triples[triple]
                pf.write(f"{triple[0]}, {triple[1]}: {sims[0]}\t")
                pf.write(f"{triple[0]}, {triple[2]}: {sims[1]}\t")
                pf.write(f"{triple[1]}, {triple[2]}: {sims[2]}\t")
                pf.write(f"{judge_side(similarity_cutoff, sims[0])}, ")
                pf.write(f"{judge_side(similarity_cutoff, sims[1])}, ")
                pf.write(f"{judge_side(similarity_cutoff, sims[2])}\n")
    
    return triples


def main(args):
    """
    Step 1: get the cutoff of the current species
    Step 2: traverse the whole file, 
            record the {<gene_segment_01, gene_segment_02>, similarity} pairs.
    Step 3: extract triple information
    """
    ### Step 1
    species_list = calculate_simi.parse_species_list(args)
    if species_list is None or len(species_list) <= 0:
        print("No valid input.")
        return

    print(f"The program will deal with the following files: {species_list}")

    # Check the output directory.
    current_dir = os.path.dirname(os.path.abspath(__file__))
    directory = "../output"
    if not os.path.exists(directory):
        os.makedirs(directory)

    # species_list = ['rainbow_trout', 'chum_salmon']
    for species in species_list:
        input_file = os.path.join(current_dir, species)
        blocks_avg_sim = calculate_simi.parse_dagchainer_output(input_file)

        similarity_cutoff = find_binormal_cutoffs(blocks_avg_sim, None, None)
        # Step 1 finished for each species

        ### Step 2
        genome_pair_sim = traverse_each_species(input_file)

        ### Step 3
        triple_output_file = os.path.join(current_dir, directory, 
                                          os.path.basename(species) + 
                                          '.triples')
        extract_triples(genome_pair_sim, similarity_cutoff, triple_output_file)


if __name__ == "__main__":
    main(sys.argv)