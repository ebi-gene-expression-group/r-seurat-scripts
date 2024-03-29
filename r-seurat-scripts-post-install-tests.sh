#!/usr/bin/env bash

script_dir=$(cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
script_name=$0

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
test_data_transfer_url='https://www.dropbox.com/s/1zxbn92y5du9pu0/pancreas_v3_files.tar.gz?dl=1'
test_working_dir=`pwd`/'post_install_tests'
test_anndata_url="https://www.dropbox.com/s/ngs3p8n2i8y33hj/pbmc3k.h5ad?dl=0"
export test_data_transfer_file=$test_working_dir/pancreas_v3_files.tar.gz
export test_data_archive=$test_working_dir/`basename $test_data_url`
export test_anndata_file=$test_working_dir/$(basename $test_anndata_url | sed 's/?dl=0//')

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
export data_dir=$test_working_dir/test_data

mkdir -p $test_working_dir
mkdir -p $output_dir
mkdir -p $data_dir

################################################################################
# Fetch test data
################################################################################

if [ ! -e "$test_data_archive" ]; then
    wget $test_data_url -P $test_working_dir
    wget $test_data_transfer_url -O $test_data_transfer_file
fi

if [ ! -e "$test_anndata_file" ]; then
    wget $test_anndata_url -O $test_anndata_file
fi

################################################################################
# List tool outputs/ inputs
################################################################################

export raw_matrix="$data_dir/matrix.mtx"
export raw_matrix_object="$output_dir/raw_matrix.rds"
export raw_seurat_object="$output_dir/raw_seurat.rds"
export raw_seurat_object_from_tab="$output_dir/raw_seurat_from_tab.rds"
export transfer_metadata_object="$data_dir/pancreas_metadata.rds"
export transfer_expression_object="$data_dir/pancreas_expression_matrix.rds"
export transfer_out_dir="$data_dir/transfer_out"
export filtered_seurat_object="$output_dir/filtered_seurat.rds"
export normalised_seurat_object="$output_dir/normalised_seurat.rds"
export variable_genes_seurat_object="$output_dir/variable_genes_seurat.rds"
export variable_genes_list="$output_dir/filtered_genes.txt"
export test_genes="$output_dir/random_genes.txt"
export scaled_seurat_object="$output_dir/scaled_seurat.rds"
export pca_seurat_object="$output_dir/pca_seurat.rds"
export pca_embeddings_file="$output_dir/pca_embeddings.csv"
export pca_loadings_file="$output_dir/pca_loadings.csv"
export pca_stdev_file="$output_dir/pca_stdev.txt"
export pca_image_file="$output_dir/pcatest.png"
export neighbours_seurat_object="$output_dir/neighbours_seurat.rds"
export cluster_seurat_object="$output_dir/cluster_seurat.rds"
export cluster_text_file="$output_dir/clusters.txt"
export tsne_seurat_object="$output_dir/tsne_seurat.rds"
export tsne_embeddings_file="$output_dir/tsne_embeddings.csv"
export html_output_dir="$output_dir/html_out"
export marker_text_file="$output_dir/markers.csv"
export anchor_object="$output_dir/anchor_object.rds"
export integrated_obj="$output_dir/integrated_object.rds"
export features_obj="$output_dir/features_object.rds"
export classify_query="$transfer_out_dir/sep_by_tech_fluidigmc1.rds"
export classify_result_object="$output_dir/classify_result.rds"
export classify_result_anchors_object="$output_dir/classify_anchorset.rds"
export pca_integrated_object="$output_dir/pca_integrated_object.rds"
export umap_result_object="$output_dir/integrated_obj_umap.rds"
export scaled_integrated_object="$output_dir/scaled_integrated_object.rds"
export pca_integrated_embeddings="$output_dir/pca_integrated_embeddings.txt"
export pca_integrated_loadings="$output_dir/pca_integrated_loadings.txt"
export pca_integrated_stdev="$output_dir/pca_integrated_stdev.txt"
export conserved_markers_result="$output_dir/conserved_markers_result.tsv"
export tmp_conserved_markers_fn_object="$output_dir/tmp_conserved_markers_fn_object.rds"
export tmp_conserved_markers_cl_object="$output_dir/tmp_conserved_markers_cl_object.rds"
export umap_map_query_result_object="$output_dir/umap_map_query_result_object.rds"
export singlecellexperiment_converted_cluster_object="$output_dir/singlecellexperiment_converted_cluster_object.rds"
export loom_converted_cluster_object="$output_dir/loom_converted_cluster_object"
export seurat_from_loom_cluster_object="$output_dir/seurat_from_loom_cluster_object.rds"
export h5seurat_cluster_object="$output_dir/h5seurat_cluster_object.h5seurat"
export h5seurat_to_seurat_cluster_object="$output_dir/h5seurat_to_seurat_cluster_object.rds"
export feature_plot_png_result="$output_dir/feature_plot_png_result.png"
export vln_plot_pdf_result="$output_dir/vln_plot_pdf_result.png"
export ridge_plot_eps_result="$output_dir/ridge_plot_eps_result.png"
export dot_plot_jpg_result="$output_dir/dot_plot_jpg_result.png"
export dim_plot_svg_result="$output_dir/dim_plot_svg_result.png"
export do_heatmap_jpg_result="$output_dir/heatmap_plot_jpj_result.jpg"
export feature_plot_rds_result="$output_dir/feature_plot_rds_result.rds"
export HoverLocator_result="$output_dir/HoverLocator_result.html"
## Test parameters- would form config file in real workflow. DO NOT use these
## as default values without being sure what they mean.

# Normalisation. See Seurat ?FilterCells

export min_genes=500
export min_umi=1000

# Normalisation. See Seurat ?NormalizeData

export normalisation_method='LogNormalize'
export scale_factor=10000

# Find variable genes. See ?FindVariableGenes

export mean_function='ExpMean'
export dispersion_function='LogVMR'
export fvg_x_low_cutoff='0.1'
export fvg_x_high_cutoff='8'
export fvg_y_low_cutoff='1'
export fvg_y_high_cutoff='Inf'

# Scale and center the data. See ?ScaleData
export vars_to_regress='nCount_RNA'
export model_use='linear'
export use_umi='TRUE'
export scale_max='10'
export block_size='1000'
export min_cells_to_block='1000'
export check_for_norm='TRUE'

# Run PCA
export pcs_compute=50
export use_imputed='FALSE'

# Plot PCA
export pca_dim_one=1
export pca_dim_two=2
export pt_size=1
export label_size=4
export weight_by_var="TRUE"
export do_label='FALSE'
export group_by='ident'
export pca_png_width=1000
export pca_png_height=1000
export pca_cols_use='red,blue,green,yellow,orange,pink,purple,black'
export pca_plot_order='7,6,5,4,3,2,1,0'

# Find clusters
export reduction_type='pca'
export dims_use='1,2,3,4,5,6,7,8,9,10'
export k_param=30
export resolution=0.8
export cluster_algorithm=1
export cluster_tmp_file_location='/tmp'

# Find neighbours
export compute_snn=TRUE

# t-SNE
export tsne_do_fast='TRUE'

# Marker detection
export logfc_threshold=0.25
export marker_min_pct=0.1
export marker_only_pos='FALSE'
export marker_test_use='wilcox'
export marker_max_cells_per_ident='Inf'
export marker_min_cells_gene=3
export marker_min_cells_group=3

# Integration
inputs_integration=""
for tech in celseq celseq2 smartseq2; do
    inputs_integration=${inputs_integration},${transfer_out_dir}/sep_by_tech_${tech}_norm_fvg.rds
done
export inputs_integration=$(echo $inputs_integration | sed s/^,// )
export integration_anchor_features=2000

# Classify against reference


################################################################################
# Test individual scripts
################################################################################

# Make the script options available to the tests so we can skip tests e.g.
# where one of a chain has completed successfullly.

export use_existing_outputs

# Derive the tests file name from the script name

tests_file="${script_name%.*}".bats

# Execute the bats tests

$tests_file
