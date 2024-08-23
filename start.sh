#!/bin/bash

# Prelude: Since the gff file is too large, compressed them
echo "*********************** Prelude: Check the gff files. ***********************\n"
if [ ! -d "data/gff" ]; then
	if [ -f "data/gff.tar.gz" ]; then
		cd data
		tar -zxvf gff.tar.gz
		cd ../
	else
		echo "No gff.tar.gz found.\n"
	fi
else
	echo "Directory data/gff exists."
fi

## Set the input folder or file.
input_path="../data/paralogs_outputs/"
# input_path="../data/paralogs outputs/Rainbow trout"
# input_path="../data/paralogs outputs/Northern pike 1"
if [ ! -z "$1" ]; then
	input_path="../$1"
else
	echo "\n****************************************************************************"
	echo "Your input path or file is empty. We will set it by default."
	echo "If you want to set it, please retry [ sh start.sh your_file_or_folder ]."
	echo "****************************************************************************\n"
fi
# echo "Finally, the input_path is set as: \n $input_path \n"

## Execute the script.
cd src
echo "\n============================= PROCESSING ============================="
echo "\n************** pair ***************"
python3 calculate_simi.py "$input_path"
echo "\n************** singleton ***************"
python3 extract_singletons.py "$input_path"
echo "\n************** triplet ***************"
python3 extract_triple_dis.py "$input_path"
cd ..

# print results
echo "\n============================= RESULTS ============================="
echo "\n************** singleton ***************"
tail output/*.singleton.t1
echo "\n************** pairs ***************"
tail output/*.similarity.t1
echo "\n************** triplets ***************"
tail output/*.triplets.t1