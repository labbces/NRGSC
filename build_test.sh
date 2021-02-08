main_dir="/home/j/Insync/jmmunozp@usp.br/Google Drive - Shared with me/NitrogenResponsiveGenotypes/Salmon2/SalmonQuant"
test_dir="/home/j/BIOINFORMATICA/test_RNAseq"

ls "$main_dir" | grep Sample > ${test_dir}/samples.txt

for i in $(cat ${test_dir}/samples.txt)
do
mkdir -p "${test_dir}"/"$i"
head -n1000 "${main_dir}"/"${i}"/quant.sf > "${test_dir}"/"${i}"/quant.sf
echo "making file Sample_"$i""
done


