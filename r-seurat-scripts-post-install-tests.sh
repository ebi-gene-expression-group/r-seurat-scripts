#!/usr/bin/env bash

# This is a test script designed to test that everything works in the various
# accessory scripts in this package. Parameters used have absolutely NO
# relation to best practice and this should not be taken as a sensible
# parameterisation for a workflow.

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
test_genes="$output_dir/random_genes.txt"
scaled_seurat_object="$output_dir/scaled_seurat.rds"

## Test parameters- would form config file in real workflow. DO NOT use these
## as default values without being sure what they mean.

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

# Scale and center the data. See ?ScaleData
vars_to_regress='nUMI'
model_use='linear'
use_umi='TRUE'
do_scale='TRUE'
do_center='TRUE'
scale_max='10'
block_size='1000'
min_cells_to_block='1000'
check_for_norm='TRUE'

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

# Get a random set of genes to use in testing argments to scale-data.R

run_command "get-random-genes.R $normalised_seurat_object $test_genes 10000" $test_genes

# Run scale-data.R

run_command "scale-data.R -i $variable_genes_seurat_object -e $test_genes -v $vars_to_regress -m $model_use -u $use_umi -s $do_scale -c $do_center -x $scale_max -b $block_size -d $min_cells_to_block -a $assay_type -n $check_for_norm -o $scaled_seurat_object" $scaled_seurat_object

################################################################################
# Finish up
################################################################################

echo "All tests passed"
exit 0