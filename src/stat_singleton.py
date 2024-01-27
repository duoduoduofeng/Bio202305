import re
import os
import sys
import json
import copy


def extract_cds(dag_chain_file):
	cds_pairs = {}

	with open(dag_chain_file, 'r') as fin:
		row_id = 0
		new_block_start = ""

		for line in fin:
			row_id += 1
			line = line.strip()

			# check data validation
			fields = line.split('\t')
			
			if line.startswith("#"):
				if len(fields) != 6:
					print(f"Invalid line {row_id} block starting.")
					continue
				
				# b28938_NC_027327.1||a28938_NC_027300.1
				pairkey = f"{fields[2]}||{fields[3]}"
				if pairkey not in cds_pairs:
					cds_pairs[pairkey] = {}
				new_block_start = line
				cds_pairs[pairkey][new_block_start] = {}
            	
			if len(fields) != 10:
				if not line.startswith("#"):
					print(f"Invalid line {row_id} with {len(fields)} fields.")
				continue

			# chr1||start1||stop1||name1||strand1||type1||db_feature_id1||genome_order1||percent_id1
			# chr2||start2||stop2||name2||strand2||type2||db_feature_id2||genome_order2||percent_id2
			chain_info = fields[1].split('||')
			chain_info_2 = fields[5].split('||')
			if len(chain_info) != 9 or len(chain_info_2) != 9:
				print(f"Invalid line {row_id} with fault \
					columns for the two genome sequence.")
				continue

			chr1 = chain_info[0]
			start1 = int(chain_info[1])
			stop1 = int(chain_info[2])

			chr2 = chain_info_2[0]
			start2 = int(chain_info_2[1])
			stop2 = int(chain_info_2[2])

			simi = float(chain_info[-1])

			chr_pair = f"{chr1}||{chr2}"
			cds_pairs[pairkey][new_block_start]["chr1"] = chr1
			cds_pairs[pairkey][new_block_start]["chr2"] = chr2
			if "cds" not in cds_pairs[pairkey][new_block_start]:
				cds_pairs[pairkey][new_block_start]["cds"] = []
			cds_pairs[pairkey][new_block_start]["cds"].append((start1, stop1, start2, stop2, simi))

	return cds_pairs


def extract_singleton(gff_file):
	gene_singletons = {}
	with open(gff_file, 'r') as fin:
		for line in fin:
			parts = line.strip().split("\t")
			if len(parts) == 9:
				chr1 = parts[0]
				typestr = parts[2]
				start = int(parts[3])
				stop = int(parts[4])

				if typestr in ["gene", "mRNA"]:
					if chr1 not in gene_singletons:
						gene_singletons[chr1] = {}
					thekey = f"{typestr}_{start}_{stop}"
					if thekey not in gene_singletons[chr1]:
						gene_singletons[chr1][thekey] = 1
					else:
						gene_singletons[chr1][thekey] += 1

	return gene_singletons


def calculate_singletons(cds_pairs, gene_singletons, simi_threshold = 87.9):
	for pair in cds_pairs:
		for chr1 in gene_singletons:
			if chr1 in pair:
				for blockkey in cds_pairs[pair]:
					if chr1 in blockkey:
						stats = {}
						# cds_pairs[pair][blockkey]["stats"] = stats # reset the statistic information
						
						if chr1 == cds_pairs[pair][blockkey]["chr1"]:
							# find the singletons on chain 1
							cds = cds_pairs[pair][blockkey]["cds"]
							if len(cds) > 1:
								for i in range(len(cds)-1):
									first_end = cds[i][1]
									second_start = cds[i+1][0]
									simi = (cds[i][-1] + cds[i+1][-1]) / 2
									simi = round(simi, 2)
									t1ort2 = "t1"
									if simi > simi_threshold:
										t1ort2 = "t2"

									gene_count = 0
									for gene in gene_singletons[chr1]:
										typestr, start, stop = gene.split("_")
										start = int(start)
										stop = int(stop)

										if start > first_end and stop < second_start:
											gene_count += 1
									if gene_count > 0:
										if "chr1" not in stats:
											stats["chr1"] = []
										stats["chr1"].append((i, i+1, simi, t1ort2, gene_count))
						
						if chr1 == cds_pairs[pair][blockkey]["chr2"]:
							# find the singletons on chain 2
							cds = cds_pairs[pair][blockkey]["cds"]
							if len(cds) > 1:
								for i in range(len(cds)-1):
									first_end = cds[i][3]
									second_start = cds[i+1][2]
									simi = (cds[i][-1] + cds[i+1][-1]) / 2
									simi = round(simi, 2)
									t1ort2 = "t1"
									if simi > simi_threshold:
										t1ort2 = "t2"

									gene_count = 0
									for gene in gene_singletons[chr1]:
										typestr, start, stop = gene.split("_")
										start = int(start)
										stop = int(stop)

										if start > first_end and stop < second_start:
											gene_count += 1
									if gene_count > 0:
										if "chr2" not in stats:
											stats["chr2"] = []
										stats["chr2"].append((i, i+1, simi, t1ort2, gene_count))

						if len(stats) > 0:
							cds_pairs[pair][blockkey]["stats"] = stats



if __name__ == "__main__":
	dag_chain_file = "../data/atlantic_salmon/atlantic_salmon_DAGChainer.txt"
	cds_pairs = extract_cds(dag_chain_file)

	# for pair in cds_pairs:
	# 	print(f"############# {pair}")

	# 	for blockkey in cds_pairs[pair]:
	# 		print(f"{blockkey}")
	# 		print(f"{json.dumps(cds_pairs[pair][blockkey])}\n")

	gff_file = "../data/atlantic_salmon/atlantic_salmon.gff"
	gene_singletons = extract_singleton(gff_file)
	# for chr1 in gene_singletons:
	# 	print(chr1)
	# 	print(f"{chr1}\t{json.dumps(gene_singletons[chr1])}")

	simi_threshold = 87.9
	calculate_singletons(cds_pairs, gene_singletons, simi_threshold)
	for pair in cds_pairs:
		print(f"\n############# {pair}")

		for blockkey in cds_pairs[pair]:
			if "stats" in cds_pairs[pair][blockkey]:
				print(f"{blockkey}")
				print(f"{json.dumps(cds_pairs[pair][blockkey])}")


