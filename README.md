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
    import pyMUZ.gen.cli as cli
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
1. Download Anaconda:
    - Mac: https://docs.anaconda.com/anaconda/install/mac-os/
    - Windows: https://docs.anaconda.com/anaconda/install/windows/
    - Linux: https://docs.anaconda.com/anaconda/install/linux/
2. Download Git: https://github.com/git-guides/install-git
3. Clone pyMUZ from github:
    ```shell
    cd ~
    mkdir git
    cd git
    git clone https://github.com/marczepeda/pyMUZ.git
    cd pyMUZ 
    ```
4. Make the environment and install pyMUZ:
    ```shell
    conda env create -f pyMUZ.yml # When conda asks you to proceed, type "y"
    conda activate pyMUZ
    pip install -e . # Include the "."
    conda deactivate
    ```
5. (Optional) Install PrimeDesign by downloading Docker: https://www.docker.com/
    ```shell
    docker pull pinellolab/primedesign # Docker desktop app needs to be open
    ```
### Update
1. Enter the environment and delete pyMUZ:
    ```shell
    cd ~
    cd git
    conda activate pyMUZ
    pip uninstall pyMUZ # Enter 'Y' when prompted
    rm -r pyMUZ # Enter 'Y' three times to completely remove the folder
    ```
2. Clone pyMUZ from github and install pyMUZ:
    ```shell
    git clone https://github.com/marczepeda/pyMUZ.git
    cd pyMUZ
    pip install -e . # Include the "."
    conda deactivate
    ```