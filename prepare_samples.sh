#!/bin/bash

# Run "source activate qiime1" prior to running this script
#
# This script takes two arguments:
#    - the project name
#    - metadata filename
# The result is a single FASTA file containing the sequences from all the
# samples listed in the metadata file, joined, quality filtered, and
# named appropriately for the following steps in the pipeline.

if ! which align_seqs.py > /dev/null 2>&1; then
    echo
    echo "Run 'source activate qiime1' first."
    echo
    exit 1;
fi

if [[ $# -ne 2 ]]; then
    echo
    echo "This script requires two arguments:"
    echo "  - the project name"
    echo "  - the filename of the metadata file"
    echo
    exit 1;
fi

project=$1
metadata_file=$2
script=`basename $0`
script_base=${script%.*}

date=$(date +%Y%m%d)

results_dir=prepared_fasta_${date}

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

# create a file to store sample statistics
stats_file=${results_dir}/scripts/${project}_stats.txt
echo "sample	16S_seqs	ITS_seqs" > ${stats_file}

# read in the sample names from the metadata file
samples=( $(awk '{print $1}' ${metadata_file}) )
samples=("${samples[@]:1}")	# remove the column header (first entry)

echo "${#samples[@]} samples listed in metadata file"
echo "samples: ${samples[@]}"
echo

samples_16S=()
samples_ITS=()
fasta_16S=()
fasta_ITS=()

function count_seqs() {
    seq_type=$1
    sample=$2
    if [[ $seq_type == "16S" ]]; then
        samples_16S+=("$sample")
        num_16S_seqs=$(grep "^>" "fasta/${sample}_${seq_type}.fasta" | wc -l)
        fasta_16S+=("${sample}/fasta/${sample}_${seq_type}.fasta")
        echo "${num_16S_seqs} sequences"
    else
        samples_ITS+=("$sample")
        num_ITS_seqs=$(grep "^>" "fasta/${sample}_${seq_type}.fasta" | wc -l)
        fasta_ITS+=("${sample}/fasta/${sample}_${seq_type}.fasta")
        echo "${num_ITS_seqs} sequences"
    fi
}

for sample in "${samples[@]}"
do
    if [[ ! -d "${sample}" ]]; then
        echo "No directory for sample ${sample}"
    else
        cd $sample
        num_16S_seqs="NA"
        num_ITS_seqs="NA"

        for seq_type in 16S ITS; do

            if [[ -d "fastq_${seq_type}" ]]; then

                if [[ -e "fasta/${sample}_${seq_type}.fasta.gz" ]]; then
                    echo "Already processed: ${sample} ${seq_type}"
                    gunzip fasta/${sample}_${seq_type}.fasta.gz
                    count_seqs ${seq_type} ${sample}
                    continue
                fi
                if [[ -e "fasta/${sample}_${seq_type}.fasta" ]]; then
                    echo "Already processed: ${sample} ${seq_type}"
                    count_seqs ${seq_type} ${sample}
                    continue
                fi

                echo "Processing sample ${sample} ${seq_type}"

                # paired end joining and quality filtering
                pipits_getreadpairslist -i fastq_${seq_type} -o list_file
                pipits_prep -i fastq_${seq_type} -o pipits_prep_output -l list_file

                if [[ ! -d "fasta" ]]; then
                    mkdir fasta
                fi

                mv pipits_prep_output/prepped.fasta fasta/${sample}_${seq_type}.fasta
                rm -rf pipits_prep_output
                rm list_file

                # add sample labels to each sequence in FASTA file
                cd fasta
                cat > mfile.txt << EOF
#SampleID	BarcodeSequence	LinkerPrimerSequence	InputFileName	Description
${sample}			${sample}_${seq_type}.fasta	none
EOF

                add_qiime_labels.py -m mfile.txt -i . -c InputFileName -o update
                mv update/combined_seqs.fna ./${sample}_${seq_type}.fasta
                rm -rf update
                rm mfile.txt
                cd ..

                count_seqs ${seq_type} ${sample}

            fi

        done

        cd ..
        echo "${sample}	${num_16S_seqs}	${num_ITS_seqs}" >> ${stats_file}
    fi
done
echo

# concatenate the individual sample FASTA files into single, large FASTA files
if [[ ${#fasta_16S[@]} -gt 0 ]]; then
    cat ${fasta_16S[@]} > ${results_dir}/${project}_16S.fasta
    nseqs=$(grep "^>" ${results_dir}/${project}_16S.fasta | wc -l)
    gzip ${results_dir}/${project}_16S.fasta
fi
echo "${#samples_16S[@]} samples with 16S data"
echo "samples: ${samples_16S[@]}"
echo "${nseqs} 16S sequences total in ${results_dir}/${project}_16S.fasta"

if [[ ${#fasta_ITS[@]} -gt 0 ]]; then
    cat ${fasta_ITS[@]} > ${results_dir}/${project}_ITS.fasta
    nseqs=$(grep "^>" ${results_dir}/${project}_ITS.fasta | wc -l)
    gzip ${results_dir}/${project}_ITS.fasta
fi
echo "${#samples_ITS[@]} samples with ITS data"
echo "samples: ${samples_ITS[@]}"
echo "${nseqs} ITS sequences total in ${results_dir}/${project}_ITS.fasta"

# go through and gzip all the FASTA files for the individual samples
for sample in "${samples[@]}"
do
    if [[ -d "${sample}/fasta" ]]; then
        cd ${sample}/fasta
        gzip *
        cd ../..
    fi
done

