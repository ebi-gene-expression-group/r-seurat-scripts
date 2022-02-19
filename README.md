# seurat-scripts 0.3.0 for Seurat 3.2.3 [![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/seurat-scripts/README.html)

In order to wrap Seurat's internal workflow in any given workflow language, it's important to have scripts to call each of those steps, which is what this package provides.

This version of seurat-scripts uses native conversions to Loom (thoroughly tested), SCE and AnnData.

## Install

The recommended method for script installation is via a Bioconda recipe called seurat-scripts.

With the [Bioconda channels](https://bioconda.github.io/#set-up-channels) configured the latest release version of the package can be installed via the regular conda install command:

```
conda install seurat-scripts
```

## Test installation

There is a test script included:

```
r-seurat-scripts-post-install-tests.sh
```

This downloads a number of datasets and executes all of the scripts described below and more.

## Commands

Currently wrapped Seurat functions are described below. Each script has usage instructions available via --help, consult function documentation in Seurat for further details.

These instructions might be currently outdated for Seurat 3.2.3, the best source of information is `--help` on each command.

###  Read10X(): read 10X data and create Seurat object from the matrix

```
seurat-read-10x.R -d <10x data directory> -o <output matrix object in .rds format>
```

### FilterCells(): filter out poor-quality cells

```
seurat-filter-cells.R -i <raw Seurat object in .rds format> -s nGene,nUMI -l <min_genes>,<min_umi> -o <output Seurat object in .rds format>
```

### NormalizeData(): normalise the expression values

```
seurat-normalise-data.R -i <filtered Seurat object in .rds format> -a <assay type> -n <normalisation method> -s <scale_factor> -o <normalised Seurat object in .rds format>
```

### FindVariableGenes(): find variable genes

```
seurat-find-variable-genes.R -i <normalised Seurat object in .rds format> -m <mean function> -d <dispersion function> -l <fvg x low cutoff> -j <fvg x high cutoff> -y <fvg y low cutoff> -z <fvg y high cutoff> -o <output Seurat object in .rds format> -t <variable genes list in text format>
```

### ScaleData(): scale expression values

```
seurat-scale-data.R -i <Seurat object with variable genes, in .rds format> -e <test genes> -v <variables to regress> -m <model to use> -u <use umi> -s <do scale> -c <do center> -x <scale max> -b <block size> -d <min cells to block> -a <assay type> -n <check for norm> -o <output Seurat object in .rds format>
```

### RunPCA(): run a principal components analysis

```
seurat-run-pca.R -i <Seurat object with scaled expression valus in .rds format> -e <test genes> -p <pcs to compute> -m <use imputed> -o <output Seurat object in .rds format> -b <pca embedding in text format> -l <pca loadings file in format> -s <pca stdev file, text format>
```

### FindTransferAnchors(): project a reference on a query (integration)

```
seurat-find-transfer-anchors.R -i <seurat object with computed dimension reduction used as query, .rds format> -r <seurat object with computed dimension reduction used as reference, .rds format> -o <anchorSet object with anchor matrix, .rds format> -n <normalization method: pcaproject or cca> -f <features to use for dimensional reduction> -d <which dimensions to use from the reduction to specify the neighbor search space, a:b format>    
```

### DimPlot(): plot dimension reductions

```
seurat-dim-plot.r -i <Seurat object with computed dimension reductions, .rds format> -r <dimension reduction, e.g. pca> -a <dim 1> -b <dim 2> -p <pt size> -l <label size> -d <do label> -f <group by> -t <title> -w <png width> -j <png height> -o <image file>
```

### RunTSNE(): run-tSNE analysis

```
seurat-run-tsne.r -i <Seurat object with computed PcA, .rds format> -r <reduction type> -d <dims to use> -e <file with genes to use> -f <do fast tsne> -o <output Seurat object in .rds format> -b <tsne embeddings in csv format>
```

### FindClusters(): generate cell clusters from expression values

```
seurat-find-clusters.r -i <<Seurat object with computed dimension reductions, .rds format>> -e <test genes> -u <dimension reduction, e.g. pca> -d <dims to use> -k <k value> -r <resolution> -a <cluster algorithm> -m <cluster tmp file location> -o <output Seurat object in .rds format> -t <clusters in txt format>
```

## Accessory scripts

### Get a random set of genes

```
seurat-get-random-genes.R <Seurat object in .rds format> <output text file> <ngenes>
```

### Convert formats

Multiple format conversion through Seurat 3:
- Possible input formats:
  - Seurat
  - AnnData (versions contemporary to Seurat 3).
  - Loom (versions contemporary to Seurat 3)
  - SingleCellExperiment
- Possible output formats:
  - Seurat
  - Loom
  - SingleCellExperiment


```
seurat-convert.R -i inputfile.rds --input-format seurat -o output.loom --output-format loom
```

This functionality should be used with care, as some elements of the objects can be lost in some conversions. When converting to and from Loom, be careful about table headers that might have offending characters for R tables.

### Export to CellBrowser

Exports a Seurat RDS object and an (optional) markers file to a format that can be read by UCSC CellBrowser:

```
seurat-export-cellbrowser.R -i <Seurat object in .rds format> [-m markers.csv] -o <directory_for_output>
```
