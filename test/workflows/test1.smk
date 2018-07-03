import os
import sys
from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider

HTTP = HTTPRemoteProvider()

test_data = 'https://s3-us-west-2.amazonaws.com/10x.files/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz'
output_dir = 'seurat_test_results'
uncompresssed_test_data_dir = '%s/test_data' % output_dir

# Main test functions

rule test_read: 
    input: "%s/seurat_raw.rds" % output_dir

# Uncompress the test data

rule uncompress_data:
    input: HTTP.remote(test_data, keep_local=True)
    output: 
        barcodes = temp('%s/barcodes.tsv' % uncompresssed_test_data_dir),  
        genes = temp('%s/genes.tsv' % uncompresssed_test_data_dir),  
        mtx = temp('%s/matrix.mtx' % uncompresssed_test_data_dir) 
    shell:
        "mkdir -p {uncompresssed_test_data_dir} && tar -xvzf {input}"
        " --strip-components 2 -C {uncompresssed_test_data_dir}" 

# Read the 10x data, testing the read-10x.R script

rule read_10x_data:
    input:
        barcodes = '%s/barcodes.tsv' % uncompresssed_test_data_dir,  
        genes = '%s/genes.tsv' % uncompresssed_test_data_dir,  
        mtx = '%s/matrix.mtx' % uncompresssed_test_data_dir 
    output:
        object = "%s/10x_data.rds" % output_dir
    params: 
        data_dir = '%s' % uncompresssed_test_data_dir
    shell:
        "read-10x.R -d {params.data_dir} -o {output.object}"  
        
# Create a Seurat object, testing create-seurat-object.R

rule create_seurat_object:
    input: 
        object = "%s/10x_data.rds" % output_dir
    output:
        object = "%s/seurat_raw.rds" % output_dir
    shell:
        "create-seurat-object.R -i {input.object} -o {output.object}"

# Cleanup

onsuccess:
    shell("rm -rf %s" % output_dir)
    print("Workflow finished, no error") 
