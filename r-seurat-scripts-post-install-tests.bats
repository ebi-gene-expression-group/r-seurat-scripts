#!/usr/bin/env bats

# Extract the test data

@test "Extract .mtx matrix from archive" {
    if [ "$use_existing_outputs" = 'true' ] && [ -f "$raw_matrix" ]; then
        skip "$raw_matrix exists and use_existing_outputs is set to 'true'"
    fi

    run rm -f $raw_matrix && tar -xvzf $test_data_archive --strip-components 2 -C $data_dir
    echo "status = ${status}"
    echo "output = ${output}"

    [ "$status" -eq 0 ]
    [ -f  "$raw_matrix" ]
}

# Create the Matrix object

@test "Matrix object creation from 10x" {
    if [ "$use_existing_outputs" = 'true' ] && [ -f "$raw_seurat_object" ]; then
        skip "$raw_matrix_object exists and use_existing_outputs is set to 'true'"
    fi

    run rm -f $raw_matrix_object && seurat-read.R -d $data_dir -o $raw_seurat_object
    echo "status = ${status}"
    echo "output = ${output}"

    [ "$status" -eq 0 ]
    [ -f  "$raw_seurat_object" ]
}

# Create mattrix object from tabular file

@test "Matrix object creation from tabular" {
    if [ "$use_existing_outputs" = 'true' ] && [ -f "$raw_seurat_object_from_tab" ]; then
        skip "$raw_matrix_object_from_tab exists and use_existing_outputs is set to 'true'"
    fi

    run rm -f $raw_matrix_object_from_tab && echo "file.copy(system.file('extdata', 'pbmc_raw.txt', package = 'Seurat'), 'test.txt')" | R --no-save && seurat-read.R -f test.txt -o $raw_seurat_object_from_tab
    echo "status = ${status}"
    echo "output = ${output}"

    [ "$status" -eq 0 ]
    [ -f  "$raw_seurat_object_from_tab" ]
}

# Run filter-cells.R

@test "Filter cells from a raw Seurat object" {
    if [ "$use_existing_outputs" = 'true' ] && [ -f "$filtered_seurat_object" ]; then
        skip "$filtered_seurat_object exists and use_existing_outputs is set to 'true'"
    fi

    run rm -f $filtered_seurat_object && seurat-filter-cells.R -i $raw_seurat_object -s nFeature_RNA,nCount_RNA -l $min_genes,$min_umi -o $filtered_seurat_object
    echo "status = ${status}"
    echo "output = ${output}"

    [ "$status" -eq 0 ]
    [ -f  "$filtered_seurat_object" ]
}

# Run normalise-date.R

@test "Normalise expression values" {
    if [ "$use_existing_outputs" = 'true' ] && [ -f "$normalised_seurat_object" ]; then
        skip "$normalised_seurat_object exists and use_existing_outputs is set to 'true'"
    fi

    run rm -f $normalised_seurat_object && seurat-normalise-data.R -i $filtered_seurat_object -n $normalisation_method -s $scale_factor -o $normalised_seurat_object
    echo "status = ${status}"
    echo "output = ${output}"

    [ "$status" -eq 0 ]
    [ -f  "$normalised_seurat_object" ]
}

# Run find-variable-genes.R

@test "Find variable genes" {
    if [ "$use_existing_outputs" = 'true' ] && [ -f "$variable_genes_list" ]; then
        skip "$variable_genes_list exists and use_existing_outputs is set to 'true'"
    fi

    run rm -f $variable_genes_list && seurat-find-variable-genes.R -i $normalised_seurat_object -m $mean_function -d $dispersion_function -l $fvg_x_low_cutoff -j $fvg_x_high_cutoff -y $fvg_y_low_cutoff -z $fvg_y_high_cutoff -o $variable_genes_seurat_object -t $variable_genes_list
    echo "status = ${status}"
    echo "output = ${output}"

    [ "$status" -eq 0 ]
    [ -f  "$variable_genes_list" ]
}

# Get a random set of genes to use in testing argments to scale-data.R

@test "Generate random gene list" {
    if [ "$use_existing_outputs" = 'true' ] && [ -f "$test_genes" ]; then
        skip "$test_genes exists and use_existing_outputs is set to 'true'"
    fi

    run rm -f $test_genes && run seurat-get-random-genes.R $normalised_seurat_object $test_genes 10000
    echo "status = ${status}"
    echo "output = ${output}"

    [ "$status" -eq 0 ]
    [ -f  "$test_genes" ]
}

# Scale expression values

@test "Scale expression values" {
    if [ "$use_existing_outputs" = 'true' ] && [ -f "$scaled_seurat_object" ]; then
        skip "$scaled_seurat_object exists and use_existing_outputs is set to 'true'"
    fi

    run rm -f $scaled_seurat_object && seurat-scale-data.R -i $variable_genes_seurat_object -e $test_genes -v $vars_to_regress -m $model_use -u $use_umi -x $scale_max -b $block_size -d $min_cells_to_block -n $check_for_norm -o $scaled_seurat_object
    echo "status = ${status}"
    echo "output = ${output}"

    [ "$status" -eq 0 ]
    [ -f  "$scaled_seurat_object" ]
}

# Run PCA

@test "Run principal component analysis" {
    if [ "$use_existing_outputs" = 'true' ] && [ -f "$pca_seurat_object" ]; then
        skip "$scaled_seurat_object exists and use_existing_outputs is set to 'true'"
    fi

    run rm -rf $pca_seurat_object && seurat-run-pca.R -i $scaled_seurat_object -e $test_genes -p $pcs_compute -o $pca_seurat_object -b $pca_embeddings_file -l $pca_loadings_file -s $pca_stdev_file
    echo "status = ${status}"
    echo "output = ${output}"

    [ "$status" -eq 0 ]
    [ -f  "$pca_seurat_object" ]
}

# Find Neighbours

@test "Run FindNeighbours" {
    if [ "$use_existing_outputs" = 'true' ] && [ -f "$neighbours_seurat_object" ]; then
        skip "$neighbours_seurat_object exists and use_existing_outputs is set to 'true'"
    fi

    run rm -rf $neighbours_seurat_object && seurat-find-neighbours.R -i $pca_seurat_object -o $neighbours_seurat_object --dims 1,2,3,4,5 --compute-snn --reduction pca
    echo "status = ${status}"
    echo "output = ${output}"

    [ "$status" -eq 0 ]
    [ -f  "$neighbours_seurat_object" ]
}

# Generate clusters

@test "Generate cell clusters from expression values" {
    if [ "$use_existing_outputs" = 'true' ] && [ -f "$cluster_text_file" ]; then
        skip "$pca_image_file exists and use_existing_outputs is set to 'true'"
    fi

    run rm -f $cluster_text_file && seurat-find-clusters.R -i $neighbours_seurat_object -r $resolution -a $cluster_algorithm -m $cluster_tmp_file_location -o $cluster_seurat_object -t $cluster_text_file
    echo "status = ${status}"
    echo "output = ${output}"

    [ "$status" -eq 0 ]
    [ -f  "$cluster_text_file" ]
}

# Run t-SNE
@test "Run-tSNE analysis" {
    if [ "$use_existing_outputs" = 'true' ] && [ -f "$tsne_seurat_object" ]; then
        skip "$tsne_seurat_object exists and use_existing_outputs is set to 'true'"
    fi

    run rm -f $tsne_seurat_object && seurat-run-tsne.R -i $pca_seurat_object -r $reduction_type -d $dims_use -e NULL -o $tsne_seurat_object -b $tsne_embeddings_file
    echo "status = ${status}"
    echo "output = ${output}"

    [ "$status" -eq 0 ]
    [ -f  "$tsne_seurat_object" ]
}

# Run t-SNE with perplexity
@test "Run-tSNE analysis with perplexity" {
    if [ "$use_existing_outputs" = 'true' ] && [ -f "$tsne_seurat_object" ]; then
        skip "$tsne_seurat_object exists and use_existing_outputs is set to 'true'"
    fi

    run rm -f $tsne_seurat_object && seurat-run-tsne.R -i $pca_seurat_object -r $reduction_type -d $dims_use -e NULL -o $tsne_seurat_object -b $tsne_embeddings_file --perplexity 20
    echo "status = ${status}"
    echo "output = ${output}"

    [ "$status" -eq 0 ]
    [ -f  "$tsne_seurat_object" ]
}

# Run marker detection

@test "Find markers for each cluster" {
    if [ "$use_existing_outputs" = 'true' ] && [ -f "$marker_text_file" ]; then
        skip "$marker_text_file exists and use_existing_outputs is set to 'true'"
    fi

    run rm -f $marker_text_file && seurat-find-markers.R -i $cluster_seurat_object -e NULL -l $logfc_threshold -m $marker_min_pct -p $marker_only_pos -t $marker_test_use -x $marker_max_cells_per_ident -c $marker_min_cells_gene -d $marker_min_cells_group -o $marker_text_file
    echo "status = ${status}"
    echo "output = ${output}"

    [ "$status" -eq 0 ]
    [ -f  "$marker_text_file" ]
}

@test "Export to CellBrowser" {
    if [ "$use_existing_outputs" = 'true' ] && [ -d "$html_output_dir" ]; then
        skip "$html_output_dir exists and use_existing_outputs is set to 'true'"
    fi

    run rm -rf $html_output_dir && seurat-export-cellbrowser.R -i $tsne_seurat_object -o $html_output_dir -n StudyTest
    echo "status = ${status}"
    echo "output = ${output}"

    [ "$status" -eq 0 ]
}

# Plot the PCA

@test "Plot dimension reduction" {
    if [ "$use_existing_outputs" = 'true' ] && [ -f "$pca_image_file" ]; then
        skip "$scaled_seurat_object exists and use_existing_outputs is set to 'true'"
    fi

    run rm -f $pca_image_file && seurat-dim-plot.R -i $cluster_seurat_object -r pca -a $pca_dim_one -b $pca_dim_two -p $pt_size -l $label_size -d $do_label -f $group_by -q $pca_plot_order -w $pca_png_width -j $pca_png_height -o $pca_image_file
    echo "status = ${status}"
    echo "output = ${output}"

    [ "$status" -eq 0 ]
    [ -f  "$pca_image_file" ]
}

# Extract the transfer test data
@test "Extract transfer test data from archive" {
    if [ "$use_existing_outputs" = 'true' ] && [ -f "$transfer_expression_object" ] && [ -f "$transfer_metadata_object" ]; then
        skip "$transfer_expression_object and $transfer_metadata_object exist and use_existing_outputs is set to 'true'"
    fi

    echo "Transfer file: $test_data_transfer_file"
    echo "Data dir: $data_dir"
    run rm -f $transfer_expression_object $transfer_metadata_object && tar -xvzf $test_data_transfer_file --strip-components 1 -C $data_dir
    echo "status = ${status}"
    echo "output = ${output}"

    [ "$status" -eq 0 ]
    [ -f "$transfer_expression_object" ]
    [ -f "$transfer_metadata_object" ]
}

@test "Split data by technology for later integration" {
    if [ "$use_existing_outputs" = 'true' ] && [ -f "$transfer_expression_split_obj" ]; then
        skip "$transfer_expression_split_obj exists and use_existing_outputs is set to true"
    fi

    echo "Transfer expression object: $transfer_expression_object"
    echo "Transfer metadata object: $transfer_metadata_object"
    run rm -f $transfer_expression_split_obj && \
        rm -rf $transfer_out_dir && mkdir -p $transfer_out_dir && \
        seurat-split-object.R -i $transfer_expression_object -o $transfer_out_dir --split-by tech -m $transfer_metadata_object

    [ "$status" -eq 0 ]
}

@test "Normalise and find variable features for integration" {
    for tech in celseq celseq2 smartseq2; do
      run seurat-normalise-data.R -i ${transfer_out_dir}/sep_by_tech_${tech}.rds -o ${transfer_out_dir}/tmp.rds &&
        seurat-find-variable-genes.R -i ${transfer_out_dir}/tmp.rds -o ${transfer_out_dir}/sep_by_tech_${tech}_norm_fvg.rds --output-text-file ${transfer_out_dir}/gene.txt &&
        rm -rf ${transfer_out_dir}/tmp.rds
    done

    [ "$status" -eq 0 ]
}

@test "Integrate data" {
    if [ "$use_existing_outputs" = 'true' ] && [ -f "$integrated_obj" ]; then
        skip "$integrated_obj exists and use_existing_outputs is set to true"
    fi

    echo $inputs_integration
    run rm -rf $integrated_obj && \
        seurat-integration.R -i $inputs_integration --anchor-features $integration_anchor_features -o $integrated_obj

    [ "$status" -eq 0 ]
}

@test "Classify against reference" {
     if [ "$use_existing_outputs" = 'true' ] && [ -f "$classify_result_object" ]; then
         skip "$classify_result_object exists and use_existing_outputs is set to true"
     fi


     run rm -rf $classify_result_object && \
         seurat-classify-against-reference.R -i $classify_query -r $integrated_obj -o $classify_result_object

     [ "$status" -eq 0 ]
}
