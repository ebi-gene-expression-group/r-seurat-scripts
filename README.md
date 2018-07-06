# Wrapper scripts for components of the Seurat toolchain

In order to wrap Seurat's internal workflow in any given workflow language, it's important to have scripts to call each of those steps. These scripts are being written here, and will improve in completeness as time progresses. 

## Install

You can just download and use the wrappers here as we develop them. But We are intending for these scripts to be available alongside the Seurat package in Bioconda. Prior to our finalising a version of this package and making it available through usual channels, it is available from our fork of the bioconda recipes using the commands below. Here we are assuming you have a healthy Bioconda install with the correct channels activate (see https://bioconda.github.io/index.html#set-up-channels). 

You may need to install conda-build:

```
conda install conda-build
```

Now you should be able to install using the following command:

```
cd <directory where you do your Git clones>
git clone git@github.com:ebi-gene-expression-group/bioconda-recipes.git
git checkout r-seurat-scripts
cd bioconda-recipes/recipes/r-seurat-scripts
conda build .
conda install --force --use-local r-seurat-scripts
```
