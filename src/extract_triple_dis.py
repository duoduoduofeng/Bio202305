import os
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
    

def main(args):
    """
    Step 1: get the cutoff
    Step 2: traverse the whole file, record the {<gene_segment_01, gene_segment_02>, similarity} pairs.
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