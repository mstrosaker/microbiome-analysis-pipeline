#!/bin/bash

# Run "source activate qiime1" prior to running this script
#
# This script takes five arguments:
#    - the project name
#    - either "16S" or "ITS" (no quotes)
#    - directory containing the output of generate_otus.sh
#    - metadata filename
#    - the maximum rarefaction depth
#    - a comma-separated list of columns in the metadata file

split_on_commas() {
    local IFS=,
    local word_list=($1)
    for word in "${word_list[@]}"; do
        echo "$word"
    done
}

if ! which align_seqs.py > /dev/null 2>&1; then
    echo
    echo "Run 'source activate qiime1' first."
    echo
    exit 1;
fi

if [[ $# -ne 6 ]]; then
    echo
    echo "This script requires five arguments:"
    echo "  - the project name"
    echo "  - the type of sequences (either 16S or ITS)"
    echo "  - the directory containing the output from generate_otus.sh"
    echo "  - the filename of the metadata file"
    echo "  - the maximum rarefaction depth (as an integer) *"
    echo "  - a comma-separated list of columns in the metadata file, for"
    echo "    comparison of alpha diversity, group significance, etc."
    echo
    echo -n " * obtain this value by running: "
    echo "'biom summarize-table -i <biom_file>'"
    echo "   and using a value somewhat smaller than the 'Min' value."
    echo
    exit 1;
fi

project=$1
otu_type=$2
otu_dir=$3
metadata_file=$4
max_rare_depth=$5
column_list=$6
script=`basename $0`
script_base=${script%.*}

date=$(date +%Y%m%d)

biom_file=${project}_${otu_type}_otu_table.biom
if [[ ${otu_type} == "ITS" ]]; then
    rep_set1="repseqs.fasta"
    rep_set2="repseqs"
else
    rep_set1="rep_set.fna"
    rep_set2="rep_set"
fi
rep_set_file=${project}_${otu_type}_${rep_set1}
results_dir=${project}_results_${otu_type}_${date}

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

# save the specified metadata file
echo "Using metadata from ${metadata_file}"
cp ${metadata_file} ${results_dir}/scripts/${script_base}_metadata.txt
metadata_file=${results_dir}/scripts/${script_base}_metadata.txt

# create a file for qiime parameters
params=${results_dir}/scripts/${script_base}_params.txt

cat > ${params} << EOF
beta_diversity:metrics bray_curtis,unweighted_unifrac,weighted_unifrac
alpha_diversity:metrics chao1,goods_coverage,observed_species,shannon,simpson,PD_whole_tree
summarize_taxa:level 2,3,4,5,6,7
EOF

# display the biom table and grab the min # of samples for use in following cmds
#cd ${otu_dir}
#echo
#biom summarize-table -i ${biom_file}
#echo
#e_value=$(biom summarize-table -i ${biom_file}  | grep "Min" | cut -d ":" -f 2)
#e_value=${e_value%.*}	# convert to int by dropping .0 at the end
#cd ..

# make pylogeny tree
echo "Making phylogeny tree"
multi_aligner=""
if [[ ${otu_type} == "ITS" ]]; then
    multi_aligner="-m muscle"
fi
align_seqs.py -i ${otu_dir}/${rep_set_file} ${multi_aligner} -o ${results_dir}/alignment
filter_alignment.py -i ${results_dir}/alignment/${project}_${otu_type}_${rep_set2}_aligned.fasta -o ${results_dir}/alignment
make_phylogeny.py -i ${results_dir}/alignment/${project}_${otu_type}_${rep_set2}_aligned_pfiltered.fasta -o ${results_dir}/${project}_${otu_type}_rep_set_tree.tre -t fasttree

phylo_tree=${results_dir}/${project}_${otu_type}_rep_set_tree.tre

# alpha diversity
echo "Calculating alpha diversity"
echo "  - alpha rarefaction"
alpha_rarefaction.py -i ${otu_dir}/${biom_file} -t ${phylo_tree} -m ${metadata_file} -o ${results_dir}/alpha_diversity -p ${params} -e ${max_rare_depth}

split_on_commas ${column_list} | while read item; do
    echo "  - compare alpha diversity: ${item}"
    compare_alpha_diversity.py -i ${results_dir}/alpha_diversity/alpha_div_collated/PD_whole_tree.txt -c ${item} -o ${results_dir}/species_significance_${item} -m ${metadata_file} -p fdr
done

# beta diversity
echo "Calculating beta diversity"
echo "  - jackknifed beta diversity"
jackknifed_beta_diversity.py -i ${otu_dir}/${biom_file} -t ${phylo_tree} -m ${metadata_file} -o ${results_dir}/jackknifed_beta_diversity -p ${params} -e ${max_rare_depth}

echo "  - beta diversity"
beta_diversity_through_plots.py -i ${otu_dir}/${biom_file} -t ${phylo_tree} -m ${metadata_file} -o ${results_dir}/beta_diversity -p ${params} -e ${max_rare_depth}

# taxonomy summaries
echo "Producing taxonomy summaries and calculating group significance"
echo "  - summarize taxa"
summarize_taxa_through_plots.py -i ${otu_dir}/${biom_file} -m ${metadata_file} -p ${params} -o ${results_dir}/taxonomy

split_on_commas ${column_list} | while read item; do
    echo "  - group significance: ${item}"

    group_significance.py -i ${otu_dir}/${biom_file} -m ${metadata_file} -o ${results_dir}/group_significance_kruskal_wallis_${item}.txt -s kruskal_wallis -c ${item}

    group_significance.py -i ${otu_dir}/${biom_file} -m ${metadata_file} -o ${results_dir}/group_significance_nonparametric_t_test_${item}.txt -s nonparametric_t_test -c ${item}
done

