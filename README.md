# Python code for Marc Zepeda
## Package Organization
- gen: input/output, data wrangling, generating plots, and statistics.
    ```shell
    import pyMUZ.gen.io as io
    import pyMUZ.gen.tidy as t
    import pyMUZ.gen.plot as p
    import pyMUZ.gen.stat as st
    import pyMUZ.gen.image as im
    import pyMUZ.gen.web as web
    ```
- bio: molecular biology & tissue culture workflows.
    ```shell
    import pyMUZ.bio.transfect as tf
    import pyMUZ.bio.ngs as ngs
    import pyMUZ.bio.fastq as f
    import pyMUZ.bio.clone as cl
    import pyMUZ.bio.pe as pe
    import pyMUZ.bio.pegLIT as pegLIT
    import pyMUZ.bio.qPCR as qPCR
    import pyMUZ.bio.genbank as gb
    ```
- dat: interacting with databases.
    ```shell
    import pyMUZ.dat.cosmic as co
    import pyMUZ.dat.cvar as cv
    import pyMUZ.dat.bls as bls
    ```

## Instructions
### Install
1. Download Anaconda.
    - Mac: https://docs.anaconda.com/anaconda/install/mac-os/
    - Windows: https://docs.anaconda.com/anaconda/install/windows/
    - Linux: https://docs.anaconda.com/anaconda/install/linux/
2. Download Git: https://github.com/git-guides/install-git
    - Check using command line terminal: git version
3. Download Docker. https://www.docker.com/
4. Make environment: write the following in a command line terminal...
    ```shell
    cd ~
    conda create --name pyMUZ_env python=3.11.5
    # When conda asks you to proceed, type "y" 
    
    conda activate pyMUZ_env
    mkdir git
    cd git
    ```
5. Download dependencies: write the following in a command line terminal...
    ```shell
    conda install pip
    conda install conda-forge::biopython
    
    pip install -U scikit-learn
    # Also, installs numpy, pandas, matplotlib, seaborn, scipy.
    
    conda install -c conda-forge statsmodels
    conda install anaconda::requests
    conda install bioconda::viennarna
    conda install conda-forge::python-levenshtein
    conda install conda-forge::adjusttext
    conda install anaconda::pillow
    conda install anaconda::git
    pip install dna-features-viewer
    ```
6. Install pyMUZ: write the following in a command line terminal...
    ```shell
    git clone https://github.com/marczepeda/pyMUZ.git
    cd pyMUZ
    pip install -e .
    # Include the "."
    
    docker pull pinellolab/primedesign
    # Docker desktop app needs to be open
    
    conda deactivate
    ```
### Update
1. Enter environment & delete pyMUZ: write the following in a command line terminal...
    ```shell
    cd ~
    cd git
    conda activate pyMUZ_env
    pip uninstall pyMUZ
    # Enter 'Y' when prompted
    
    rm -r pyMUZ
    # Enter 'Y' three times to completely remove the folder
    ```
2. Install pyMUZ: write the following in a command line terminal...
    ```shell
    git clone https://github.com/marczepeda/pyMUZ.git
    cd pyMUZ
    pip install -e .
    # Include the "."

    conda deactivate
    ```