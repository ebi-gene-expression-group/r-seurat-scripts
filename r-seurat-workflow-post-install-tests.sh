#!/usr/bin/env bash

test_data_url='https://s3-us-west-2.amazonaws.com/10x.files/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz'
test_data_archive=`basename $test_data_url`
test_working_dir='post_install_tests'

mkdir -p $test_working_dir
cd $test_working_dir
mkdir -p $test_working_dir/outputs

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

################################################################################
# Test individual scripts
################################################################################

# Run read-10x.R

read-10x.R -d test_data -o $raw_matrix_object
report_status read-10x.R $?

# Run create-seurat-object.R

create-seurat-object.R -i $raw_matrix_object -o $raw_seurat_object
report_status create-seurat-object.R $?

################################################################################
# Finish up
################################################################################

echo "All tests passed"
echo "Cleaning up..."
cd ..
rm -rf $test_working_dir
exit 0
