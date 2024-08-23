import sys
import os
import calculate_simi

def parse_and_extract(file_path, out_file_path):
    with open(file_path, 'r') as file, \
        open(out_file_path, 'w') as fout:
        # global initialization
        row_id = 0
        genome_default = "CURR"
        genome_default_2 = "CURR_2"

        # block initialization
        current_genome_seq = genome_default
        current_genome_seq_2 = genome_default_2
        block_count = 0

        fout.write(f"# genome_seq\tstart_1\tend_1\ttype_1\tkind_1\tfid_1\tpid_1\tgenome_seq_2\tstart_2\tend_2\ttype_2\tkind_2\tfid_2\tpid_2\n")

        for line in file:
            row_id += 1
            line = line.strip()
            
            if not line:
                continue

            # If starts with "#", new block starts.
            # Thus calculate the average of last block.
            if line.startswith("#"):
                # Record the second line
                if row_id == 2:
                    fout.write(f"{line}\t{-1}\n")
                
                if block_count <= 0:
                    # print(f"Skip line {row_id} which should be the starting.")
                    continue
                
                # record block start line and last block's size for checking
                fout.write(f"{line}\t{block_count}\n")
                    
                # re-init
                current_genome_seq = genome_default
                current_genome_seq_2 = genome_default_2
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
            start_1 = chain_info[1]
            end_1 = chain_info[2]
            type_1 = chain_info[4]
            kind_1 = chain_info[5]
            fid_1 = chain_info[6]
            pid_1 = chain_info[8]

            current_genome_seq_2 = genome_seq_2
            start_2 = chain_info_2[1]
            end_2 = chain_info_2[2]
            type_2 = chain_info_2[4]
            kind_2 = chain_info_2[5]
            fid_2 = chain_info_2[6]
            pid_2 = chain_info_2[8]

            block_count += 1

            fout.write(f"{current_genome_seq}\t{start_1}\t{end_1}\t{type_1}\t{kind_1}\t{fid_1}\t{pid_1}\t{current_genome_seq_2}\t{start_2}\t{end_2}\t{type_2}\t{kind_2}\t{fid_2}\t{pid_2}\n")


def parse_gff(file_path, out_file_path):
    with open(file_path, 'r') as file, \
        open(out_file_path, 'w') as fout:
        row_id = 0

        # head of the file
        fout.write(f"# genome_name\tseq_type\tstart\tend\treverse\tfid\n")

        for line in file:
            row_id += 1
            line = line.strip()
            
            if not line:
                continue

            # Ignore lines that start with "#"
            if line.startswith("#"):
                continue

            # "genome_name\tseq_type\tstart\tend\treverse\tfid\n"
            fields = line.split("\t")
            genome_name = fields[0]
            seq_type = fields[2]

            # chromosome has no meaning here.
            if seq_type in ["chromosome", "contig"]:
                continue

            start = fields[3]
            end = fields[4]
            reverse = 1 if fields[6] == "+" else -1

            subfields = fields[-1].split(";")
            fid = "";
            for i in range(len(subfields) - 1, 0, -1):
                if "coge_fid=" in subfields[i]:
                    fid = subfields[i].split("=")[-1]
                    break
            
            fout.write(f"{genome_name}\t{seq_type}\t{start}\t{end}\t{reverse}\t{fid}\n")


def extract(dag_file, gff_file, inpair_gff_out, notinpair_gff_out):
    with open(dag_file, 'r') as df, \
            open(gff_file, 'r') as gf, \
            open(inpair_gff_out, 'w') as inout, \
            open(notinpair_gff_out, 'w') as nout:
        
        # Preparing dag dict.
        dag_dict = {}
        row_id = 0
        for line in df:
            row_id += 1
            line = line.strip()
            
            if not line or line.startswith("#"):
                continue

            fields = line.split("\t")
            if len(fields) != 14:
                # print(f"DAG, Skip line {row_id} which doesn't contain 14 fields.")
                continue

            # sequence 1
            genome_seq = fields[0]
            if genome_seq not in dag_dict:
                dag_dict[genome_seq] = {}
            fid = fields[5]
            if fid in dag_dict[genome_seq]:
                # print(f"DAG, Skip line {row_id} which duplicated fid {fid}.")
                continue
            start = int(fields[1])
            end = int(fields[2])
            pid = float(fields[6])
            dag_dict[genome_seq][fid] = [start, end, pid]

            # sequence 2
            genome_seq = fields[7]
            if genome_seq not in dag_dict:
                dag_dict[genome_seq] = {}
            fid = fields[12]
            if fid in dag_dict[genome_seq]:
                # print(f"DAG, Skip line {row_id} which duplicated fid {fid} on the second genome sequence.")
                continue
            start = int(fields[8])
            end = int(fields[9])
            pid = float(fields[13])
            dag_dict[genome_seq][fid] = [start, end, pid]

        # for seq in dag_dict:
        #     print(f"{seq}\n{dag_dict[seq]}\n")

        # Start extracting
        row_id = 0
        inout.write(f"# genome_name\tseq_type\tstart\tend\treverse\tfid\n")
        nout.write(f"# genome_name\tseq_type\tstart\tend\treverse\tfid\n")
        for line in gf:
            row_id += 1
            line = line.strip()
            
            if not line or line.startswith("#"):
                continue

            fields = line.split("\t")
            if len(fields) != 6:
                # print(f"GFF file, skip line {row_id} which doesn't contain 6 fields.")
                continue
            
            # genome_name   seq_type        start   end     reverse fid
            genome_name = fields[0]
            if genome_name not in dag_dict:
                # print(f"GFF file, skip line {row_id} whose genome_name is not in dag dict.")
                continue
            fid_dict = dag_dict[genome_name]
            fid = fields[5]
            # inside a pair
            if fid in fid_dict:
                inout.write(f"{line}\n")
            else:
                nout.write(f"{line}\n")


def stat_singleton(dag_file, gff_file, notinpair_gff_out):
    with open(dag_file, 'r') as df, \
            open(gff_file, 'r') as gf, \
            open(notinpair_gff_out, 'w') as nout:
        
        # Preparing dag dict.
        dag_dict = {}
        row_id = 0
        blockid = 0
        for line in df:
            row_id += 1
            line = line.strip()
            
            if not line:
                continue
            if line.startswith("#"):
                blockid += 1
                continue

            fields = line.split("\t")
            if len(fields) != 14:
                # print(f"Skip line {row_id} which doesn't contain 14 fields.")
                continue

            # sequence 1
            genome_seq = fields[0]
            if genome_seq not in dag_dict:
                dag_dict[genome_seq] = {}
            fid = fields[5]
            if fid in dag_dict[genome_seq]:
                # print(f"Skip line {row_id} which duplicated fid {fid}.")
                continue
            start = int(fields[1])
            end = int(fields[2])
            pid = float(fields[6])
            dag_dict[genome_seq][fid] = [start, end, pid, blockid]

            # sequence 2
            # genome_seq = fields[7]
            # if genome_seq not in dag_dict:
            #     dag_dict[genome_seq] = {}
            # fid = fields[13]
            # if fid in dag_dict[genome_seq]:
            #     print(f"Skip line {row_id} which duplicated fid {fid} on the second genome sequence.")
            #     continue
            # start = int(fields[8])
            # end = int(fields[9])
            # pid = float(fields[13])
            # dag_dict[genome_seq][fid] = [start, end, pid]

        # Start extracting
        last_paired_fid_dict = {}
        thefid = "LAST"
        lastfid = "LAST"
        last_batch = []
        current_batch = []
        row_id = 0
        
        for line in gf:
            row_id += 1
            line = line.strip()
            
            if not line or line.startswith("#"):
                continue

            fields = line.split("\t")
            if len(fields) != 6:
                # print(f"GFF file, skip line {row_id} which doesn't contain 6 fields.")
                continue
            
            # genome_name   seq_type        start   end     reverse fid
            genome_name = fields[0]
            if genome_name not in dag_dict:
                # print(f"GFF file, skip line {row_id} whose genome_name is not in dag dict.")
                continue
            fid_dict = dag_dict[genome_name]
            fid = fields[5]
            
            # inside a pair
            if fid in fid_dict:
                last_paired_fid_dict[genome_name] = fid
                # last_batch = current_batch
                # current_batch = []
            else:
                if genome_name in last_paired_fid_dict:
                    thefid = last_paired_fid_dict[genome_name]
                    pid = dag_dict[genome_name][thefid][2]
                    blockid = dag_dict[genome_name][thefid][3]
                else:
                    pid = "None"
                    blockid = "None"
                # nout.write(f"{line}\t{thefid}\t{pid}\n")
                
                if lastfid != thefid:
                    # print(f"{lastfid}\t{thefid}")
                    # if len(last_batch) > 0:
                    lastpid = "None"
                    lastblockid = "None"
                    if lastfid in dag_dict[genome_name]:
                        lastpid = dag_dict[genome_name][lastfid][2]
                        lastblockid = dag_dict[genome_name][lastfid][3]
                    
                    for theline in last_batch:
                        nout.write(f"{theline}\t{lastfid}\t{lastpid}\t{lastblockid}\n")
                    last_batch = current_batch
                    current_batch = []
                    lastfid = thefid
                else:
                    current_batch.append(f"{line}\t{thefid}\t{pid}\t{blockid}")
        
        # clear the last batch
        for theline in current_batch:
            nout.write(f"{theline}\tNone\tNone\tNone\n")


def calculate_t(notinpair_gff_out, singleton_stat_file, singletons_between_file, cutoff):
    # not_same_block_count = 0
    singleton_between_count = 0
    t1_count = 0
    t2_count = 0

    with open(notinpair_gff_out, 'r') as fin, \
        open(singleton_stat_file, 'w') as fout, \
        open(singletons_between_file, 'w') as fout2: 
        row_id = 0
        for line in fin:
            row_id += 1
            line = line.strip()
            
            if not line:
                continue

            fields = line.split("\t")
            if len(fields) != 12:
                # print(f"In notinpair_gff_out, skip line {row_id} which doesn't contain 12 fields: {line}")
                continue

            # Filter singleton between blocks
            blockid1 = fields[8]
            blockid2 = fields[11]
            # if blockid1 != "None" and blockid2 != "None" and blockid1 != blockid2:
                # not_same_block_count += 1
            if blockid1 == "None" or blockid2 == "None" or blockid1 != blockid2:
                singleton_between_count += 1
                fout2.write(f"{line}\n")
                continue

            pid1 = fields[7]
            pid2 = fields[10]
            # if pid1 == "None":
            #     pid1 = pid2
            # if pid2 == "None":
            #     pid2 = pid1
                
            pid1 = float(pid1)
            pid2 = float(pid2)

            avg_pid = (pid1 + pid2) / 2.0
            if avg_pid <= cutoff:
                t1_count += 1
            else:
                t2_count += 1
        
        # fout.write(f"Singleton: not_same_block_count = {not_same_block_count}, \n")
        fout.write(f"Singleton: singleton_between_count = {singleton_between_count}, \n")
        fout.write(f"t1_count = {t1_count}, t2_count = {t2_count}\n")


def parse_cutoff(cutoff_para_file):
    cutoff = -1.0
    with open(cutoff_para_file, 'r') as fin:
        for line in fin:
            if line.strip().startswith("Best cutoff: "):
                cutoff = float(line.strip().split(": ")[-1])
                break
    return cutoff



def executes(files):
    # parse dag file
    dag_input_file = files[0]
    dag_output_file = files[1]
    parse_and_extract(dag_input_file, dag_output_file)

    # parse gff file
    gff_input_file = files[2]
    gff_output_file = files[3]
    parse_gff(gff_input_file, gff_output_file)

    # extract
    inpair_gff_out = f"{gff_output_file}.in"
    notinpair_gff_out = f"{gff_output_file}.not"
    extract(dag_output_file, gff_output_file, inpair_gff_out, notinpair_gff_out)

    # for statistic
    notinpair_gff_out = f"{gff_output_file}.forstat"
    stat_singleton(dag_output_file, gff_output_file, notinpair_gff_out)

    cutoff_para_file = files[4]
    cutoff = parse_cutoff(cutoff_para_file)
    if cutoff == -1:
        print(f"No valid cutoff input from file {cutoff_para_file}, skip.")
        return
    # cutoff = 86.6
    singleton_stat_file = files[5]
    singletons_between_file = files[6]
    calculate_t(notinpair_gff_out, singleton_stat_file, singletons_between_file, cutoff)


def main(args):
    species_list = calculate_simi.parse_species_list(args)
    if species_list is None or len(species_list) <= 0:
        print("No valid input.")
        return

    print(f"\n===================  Extract singletons.  ===================")
    print(f"The program will deal with the following files: \n{species_list}.\n")

    # Check the output directory.
    current_dir = os.path.dirname(os.path.abspath(__file__))
    directory = "../output"
    if not os.path.exists(directory):
        os.makedirs(directory)

    # species_list = ['rainbow_trout', 'chum_salmon']
    for species in species_list:
        dag_input_file = os.path.join(current_dir, species)
        gff_input_file = os.path.join(current_dir, "../data/gff/", 
                                      os.path.basename(species) + ".gff")
        
        print(f"****** Start dealing with {os.path.basename(species)}. ******")
        if not os.path.isfile(gff_input_file):
            print(f"No gff file for {species}, skip.")
            continue
        
        dag_output_file = os.path.join(current_dir, directory, 
                                       os.path.basename(species) + 
                                       '.dag')
        gff_output_file = os.path.join(current_dir, directory, 
                                       os.path.basename(species) + 
                                       '.gff')
        binormal_cutoff_paras_files = os.path.join(current_dir, directory, 
                                                   os.path.basename(species) + 
                                                   '.binormal_cutoff.parameters')
        singleton_stat_file = os.path.join(current_dir, directory, 
                                            os.path.basename(species) + 
                                            '.singleton.t1')
        singletons_between_file = os.path.join(current_dir, directory, 
                                            os.path.basename(species) + 
                                            '.singleton_between')
        files = [dag_input_file, dag_output_file, gff_input_file, 
                 gff_output_file, binormal_cutoff_paras_files, 
                 singleton_stat_file, singletons_between_file]
        executes(files)


if __name__ == "__main__":
    # python3 extract_sigletons.py ../data/paralogs\ outputs/Myxocyprinus\ asiaticus_ tmp.dag ../data/gff/Myxocyprinus_asiaticus.gff tmp.gff
    main(sys.argv)