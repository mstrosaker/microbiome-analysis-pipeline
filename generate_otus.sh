#!/bin/bash

# Run "source activate qiime1" prior to running this script
#
# This script takes four arguments:
#    - the project name
#    - "16S" or "ITS"
#    - "openref" or "closedref"
#    - the directory containing the output of prepare_samples.sh
#
# The results are in the OTUs_<16S|ITS>_<open|closed>ref_<date> directory

if ! which align_seqs.py > /dev/null 2>&1; then
    echo
    echo "Run 'source activate qiime1' first."
    echo
    exit 1;
fi

project=$1
seq_type=$2
ref=$3
indir=$4
script=`basename $0`
script_base=${script%.*}

if [[ $# -ne 4 ]]; then
    echo
    echo "This script requires four arguments:"
    echo "  - the project name"
    echo "  - \"16S\" or \"ITS\""
    echo "  - \"openref\" or \"closedref\" (must be closedref for ITS sequences)"
    echo "  - the directory containing the output of prepare_samples.sh"
    echo
    exit 1;
fi

if [[ $seq_type != "ITS" && $seq_type != "16S" ]]; then
    echo
    echo "The second argument must be \"16S\" or \"ITS\"."
    echo
fi

if [[ $seq_type == "ITS" && $ref != "closedref" ]]; then
    echo
    echo "The third argument must be \"closedref\" if the second is \"ITS\"."
    echo
fi

if [[ $ref != "closedref" && $ref != "openref" ]]; then
    echo
    echo "The third argument must be \"closedref\" or \"openref\"."
    echo
    exit 1;
fi

if [[ ! -d $indir ]]; then
    echo
    echo "The fourth argument must specify a directory with the output from"
    echo "prepare_samples.sh."
    echo
    exit 1;
fi

date=$(date +%Y%m%d)
results_dir=${project}_OTUs_${seq_type}_${ref}_${date}

: "${GG_OTU_BASE:?The GG_OTU_BASE environment variable should be set to the base GreenGenes database path}"
OTU_PERCENT=97
GG_OTU_REP=${GG_OTU_BASE}/rep_set/${OTU_PERCENT}_otus.fasta
GG_OTU_TAX=${GG_OTU_BASE}/taxonomy/${OTU_PERCENT}_otu_taxonomy.txt

if [[ -e ${results_dir} ]]; then
    echo
    echo "Output directory (${results_dir}) already exists."
    echo
    exit 1;
fi
echo "Results will be written to ${results_dir}"

# create results directory and copy this script to it, in case we need to
# refer to it later
mkdir ${results_dir}
mkdir ${results_dir}/scripts
cp "$(readlink -f $0)" "${results_dir}/scripts"

# save the arguments that this was invoked with to the results directory
this_pid=${BASHPID}
cat /proc/${this_pid}/cmdline | tr "\0" " " > ${results_dir}/scripts/${script_base}_cmdline.txt

# output to both the console and a log file in the results directory
log_file=${results_dir}/scripts/${script_base}.log
echo "Logging to ${log_file}"
exec > >(tee -a ${log_file} )
exec 2> >(tee -a ${log_file} >&2)

if [[ $seq_type == "16S" ]]; then
    # write a parameters file to the results directory
    params=${results_dir}/scripts/${script_base}_params.txt

    cat > params.txt << EOF
pick_otus:enable_rev_strand_match True
EOF
fi

# create a file to store statistics
stats_file=${results_dir}/scripts/${project}_${seq_type}_pick_otus_${ref}_stats.txt

if [[ $seq_type == "16S" ]]; then
    infile=${indir}/${project}_16S_chimeras_filtered.fasta
    #infile=${indir}/${project}_16S.fasta
    mkdir temp
else
    infile=${indir}/${project}_ITS.fasta
fi

gunzip ${infile}.gz
cp ${infile} ${indir}/${project}_16S.fasta
gzip ${infile}
infile=${indir}/${project}_16S.fasta

n_seqs=$(grep ">" $infile | wc -l)
echo "$n_seqs input sequences"

if [[ $seq_type == "16S" && $ref == "closedref" ]]; then

    echo "Picking OTUs"
    pick_otus.py -i ${infile} -o temp/pick_otus -r ${GG_OTU_REP} -m usearch61_ref --suppress_new_clusters --enable_rev_strand_match

    N_OTUS=$(cat temp/pick_otus/${project}_${seq_type}_otus.txt | wc -l)
    echo "    - $N_OTUS OTUs identified"

    cp temp/pick_otus/${project}_${seq_type}_otus.txt ${results_dir}/${project}_16S_otu_map.txt

    echo "Picking representative set"
    pick_rep_set.py -i temp/pick_otus/${project}_${seq_type}_otus.txt -o ${results_dir}/${project}_16S_rep_set.fna -f ${infile}

    echo "Creating OTU abundance table"
    make_otu_table.py -i temp/pick_otus/${project}_${seq_type}_otus.txt -t ${GG_OTU_TAX} -o ${results_dir}/${project}_16S_otu_table.biom

    echo "Assigning taxonomy"
    assign_taxonomy.py -o temp/taxonomy -i ${results_dir}/${project}_16S_rep_set.fna

    biom add-metadata -i ${results_dir}/${project}_16S_otu_table.biom --observation-metadata-fp temp/taxonomy/${project}_16S_rep_set_tax_assignments.txt -o ${results_dir}/${project}_16S_otu_table.biom --sc-separated taxonomy --observation-header OTUID,taxonomy

    rm -rf temp

elif [[ $seq_type == "16S" && $ref == "openref" ]]; then

    echo "Picking OTUs"
    pick_open_reference_otus.py -i ${infile} -o temp/pick_otus -r ${GG_OTU_REP} -p ${params}

    cp temp/pick_otus/final_otu_map_mc2.txt ${results_dir}/${project}_16S_openref_otu_map.txt
    cp temp/pick_otus/otu_table_mc2_w_tax.biom ${results_dir}/${project}_16S_openref_otu_table.biom
    cp temp/pick_otus/rep_set.fna ${results_dir}/${project}_16S_openref_rep_set.fna
    cp temp/pick_otus/new_refseqs.fna ${results_dir}/${project}_16S_openref_new_refseqs.fna
    cp temp/pick_otus/rep_set.tre ${results_dir}/${project}_16S_openref_rep_set.tre

    rm -rf temp

else # ITS

    echo "Extracting fungal ITS regions from the reads"
    pipits_funits -i ${infile} -o pipits_funits_output -x ITS2

    echo "Creating OTU abundance tables and assigning taxonomy"
    pipits_process -i pipits_funits_output/ITS.fasta -o ${results_dir} --Xmx 4G
    #rm -rf pipits_funits_output

    for fname in ${results_dir}/*.*
    do
        fname_new=`echo $fname | sed "s/${results_dir}\//${results_dir}\/${project}_ITS_/"`
        mv $fname $fname_new
    done

    cp -r pipits_funits_output ${results_dir}

fi

# write out the statistics file
biom summarize-table -i ${results_dir}/${project}_16S_otu_table.biom >> ${stats_file}

echo
biom summarize-table -i ${results_dir}/${project}_16S_otu_table.biom
echo

