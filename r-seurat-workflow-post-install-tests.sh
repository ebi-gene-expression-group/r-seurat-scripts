#!/usr/bin/env bash

function usage {
    echo "usage: r-seurat-workflow-post-install-tests.sh [action]"
    echo "  -action: what action to take, 'test' or 'clean'"
    exit 1
}

if [ $# -eq 0 ]; then
    echo "No arguments supplied"
    usage
fi

action=$1

test_data_url='https://s3-us-west-2.amazonaws.com/10x.files/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz'
test_data_archive=`basename $test_data_url`
test_working_dir='post_install_tests'

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

mkdir -p $test_working_dir
mkdir -p $test_working_dir/outputs

cd $test_working_dir

################################################################################
# Fetch test data 
################################################################################

if [ ! -e "$test_data_archive" ]; then
    wget $test_data_url
fi

echo "Extracting test data from archive"
mkdir -p test_data && tar -xvzf $test_data_archive --strip-components 2 -C test_data

################################################################################
# Accessory functions
################################################################################

function report_status() {
    script=$1
    status=$2

    if [ $status != 0 ]; then
        echo "FAIL: $script"
        exit 1
    else
        echo "SUCCESS: $script"
    fi
}

################################################################################
# List tool outputs/ inputs
################################################################################

raw_matrix_object='outputs/raw_matrix.rds'
raw_seurat_object='outputs/raw_seurat.rds'
filtered_seurat_object='outputs/filtered_seurat.rds'

################################################################################
# Test individual scripts
################################################################################

# Run read-10x.R

read-10x.R -d test_data -o $raw_matrix_object
report_status read-10x.R $?

# Run create-seurat-object.R

create-seurat-object.R -i $raw_matrix_object -o $raw_seurat_object
report_status create-seurat-object.R $?

# Run filter-cells.R

filter-cells.R -i $raw_seurat_object -s nGene,nUMI -l 500,1000 -o $filtered_seurat_object
report_status filter-cells.R

################################################################################
# Finish up
################################################################################

echo "All tests passed"
exit 0
