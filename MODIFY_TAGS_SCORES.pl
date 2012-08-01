
# to get only first, second and third column of bed file:
awk '{print $1"\t"$2"\t"$3}' dummy.bed > dummy_2.bed


#scores (fifth column) in tags  file_A.bed or file_B.bed are not integers where they must be. so I cannot convert them to bigbed files. I use the following awk command to transfom the scores into integers ranging from zero to 1000 and then the oupt file from this command can be converted to bigbed file.
#Note that for conversion at the moment I have only multiplied by 100 and then round it to nearest 10 but you may need to do a better conversion.

awk '{$5=-100*$5;$5 = sprintf("%.0f", $5);print $0}' file_A.bed > dummy_file_A.bed


