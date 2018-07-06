#!/usr/bin/env bash

function usage {
    echo "usage: r-seurat-workflow-post-install-tests.sh [action] [use_existing_outputs]"
    echo "  - action: what action to take, 'test' or 'clean'"
    echo "  - use_existing_outputs, 'true' or 'false'"
    exit 1
}

action=${1:-'test'}
use_existing_outputs=${2:-'false'}

if [ "$action" != 'test' ] && [ "$action" != 'clean' ]; then
    echo "Invalid action"
    usage
fi

if [ "$use_existing_outputs" != 'true' ] && [ "$use_existing_outputs" != 'false' ]; then
    echo "Invalid value ($use_existing_outputs) for 'use_existing_outputs'"
    usage
fi

test_data_url='https://s3-us-west-2.amazonaws.com/10x.files/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz'
test_working_dir=`pwd`/'post_install_tests'
test_data_archive=$test_working_dir/`basename $test_data_url`

# Clean up if specified

if [ "$action" = 'clean' ]; then
    echo "Cleaning up $test_working_dir ..."
    rm -rf $test_working_dir
    exit 0
elif [ "$action" != 'test' ]; then
    echo "Invalid action '$action' supplied"
    exit 1
fi 

# Initialise directories

output_dir=$test_working_dir/outputs
data_dir=$test_working_dir/test_data

mkdir -p $test_working_dir
mkdir -p $output_dir
mkdir -p $data_dir

cd $test_working_dir

################################################################################
# Fetch test data 
################################################################################

if [ ! -e "$test_data_archive" ]; then
    wget $test_data_url
fi

################################################################################
# Accessory functions
################################################################################

function report_status() {
    script=$1
    status=$2

    if [ $status -ne 0 ]; then
        echo "FAIL: $script"
        exit 1
    else
        echo "SUCCESS: $script"
    fi
}

# Run a command, checking the primary output depending on the value of
# 'use_existing_outputs'

run_command() {
    command=$1
    test_output=$2

    echo "$command"
    command_name=`echo "$command" | awk '{print $1}'`

    if [ -e "$test_output" ] && [ "$use_existing_outputs" == "true" ]; then
        echo "Using cached output for $command_name"
    else
        eval $command
        report_status $command_name $?
    fi
}

################################################################################
# List tool outputs/ inputs
################################################################################

raw_matrix='test_data/matrix.mtx'
raw_matrix_object="$output_dir/raw_matrix.rds"
raw_seurat_object="$output_dir/raw_seurat.rds"
filtered_seurat_object="$output_dir/filtered_seurat.rds"
normalised_seurat_object="$output_dir/normalised_seurat.rds"
variable_genes_seurat_object="$output_dir/variable_genes_seurat.rds"
variable_genes_list="$output_dir/filtered_genes.txt"

## Test parameters- would form config file in real workflow

# Normalisation. See Seurat ?FilterCells

min_genes=500
min_umi=1000

# Normalisation. See Seurat ?NormalizeData

assay_type='RNA'
normalisation_method='LogNormalize'
scale_factor=10000

# Find variable genes. See ?FindVariableGenes

mean_function='ExpMean'
dispersion_function='LogVMR'
fvg_x_low_cutoff='0.1'
fvg_x_high_cutoff='8'
fvg_y_low_cutoff='1'
fvg_y_high_cutoff='Inf'

################################################################################
# Test individual scripts
################################################################################

# Extract the test data

echo "Extracting test data from archive"
run_command "tar -xvzf $test_data_archive --strip-components 2 -C test_data" $raw_matrix

# Run read-10x.R

run_command "read-10x.R -d test_data -o $raw_matrix_object" $raw_matrix_object

# Run create-seurat-object.R

run_command "create-seurat-object.R -i $raw_matrix_object -o $raw_seurat_object" $raw_seurat_object

# Run filter-cells.R

run_command "filter-cells.R -i $raw_seurat_object -s nGene,nUMI -l $min_genes,$min_umi -o $filtered_seurat_object" $filtered_seurat_object

# Run normalise-date.R

run_command "normalise-data.R -i $filtered_seurat_object -a $assay_type -n $normalisation_method -s $scale_factor -o $normalised_seurat_object" $normalised_seurat_object

# Run find-variable-genes.R

run_command "find-variable-genes.R -i $normalised_seurat_object -m $mean_function -d $dispersion_function -l $fvg_x_low_cutoff -j $fvg_x_high_cutoff -y $fvg_y_low_cutoff -z $fvg_y_high_cutoff -o $variable_genes_seurat_object -t $variable_genes_list" $variable_genes_list

################################################################################
# Finish up
################################################################################

echo "All tests passed"
exit 0
