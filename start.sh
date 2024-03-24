#!/bin/bash

## Set the input folder or file.
# input_path="../data/paralogs outputs/"
# input_path="../data/paralogs outputs/Rainbow trout"
input_path="../data/paralogs outputs/Northern pike 1"
if [ ! -z "$1" ]; then
	input_path="../$1"
else
	echo "Your input path or file is empty. We will set it by default."
	echo "\n***************************************************"
	echo "If you want to set it, please retry [ sh start.sh your_file_or_folder ]."
	echo "***************************************************\n"
fi
echo "Finally, the input_path is set as: \n $input_path \n"

## Execute the script.
cd src
python3 calculate_simi.py "$input_path"
cd ..
