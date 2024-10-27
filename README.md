# Python code for Marc Zepeda
## Package Organization
- gen: input/output, data wrangling, generating plots, statistics, and interacting with websites.
    - import pyMUZ.gen.io as io
    - import pyMUZ.gen.tidy as t
    - import pyMUZ.gen.plot as p
    - import pyMUZ.gen.stat as st
    - import pyMUZ.gen.web as web
- bio: molecular biology workflows including cloning and sequencing.
    - import pyMUZ.bio.fastq as f
    - import pyMUZ.bio.clone as cl
    - import pyMUZ.bio.pe as pe
    - import pyMUZ.bio.pegLIT as pegLIT
    - import pyMUZ.bio.qPCR as qPCR
- dat: interacting with database APIs.
    - import pyMUZ.dat.cosmic as co
    - import pyMUZ.dat.cvar as cv
    - import pyMUZ.dat.bls as bls

## Instructions
### Clone (Option 1)
1. Download Anaconda.
    - Mac: https://docs.anaconda.com/anaconda/install/mac-os/
    - Windows: https://docs.anaconda.com/anaconda/install/windows/
    - Linux: https://docs.anaconda.com/anaconda/install/linux/
2. Download Git: https://github.com/git-guides/install-git
    - Check using command line terminal: git version
3. Write the following in a command line terminal:
    - cd ~
    - conda create --name emds_env python=3.11.5
        - When conda asks you to proceed, type "y" 
    - conda activate pyMUZ_env
    - mkdir git_pyMUZ
    - cd git_pyMUZ
    - git clone https://github.com/marczepeda/pyMUZ.git
    - pip install -e .
4. Optional: download Docker at https://www.docker.com/
    - docker pull pinellolab/primedesign
        - Docker desktop app needs to be open

### Install (Option 2)
1. Download Anaconda.
    - Mac: https://docs.anaconda.com/anaconda/install/mac-os/
    - Windows: https://docs.anaconda.com/anaconda/install/windows/
    - Linux: https://docs.anaconda.com/anaconda/install/linux/
2. Download Git: https://github.com/git-guides/install-git
    - Check using command line terminal: git version
3. Write the following in a command line terminal:
    - conda create --name pyMUZ_env python=3.11.5
        - When conda asks you to proceed, type "y" 
    - conda activate pyMUZ_env
    - pip install git+https://github.com/marczepeda/pyMUZ.git
        - Sign into github when asked.
        - Wait at least a minute for the authentication to process.
    - conda list
        - check for pyMUZ package
4. Optional: download Docker at https://www.docker.com/
    - docker pull pinellolab/primedesign
        - Docker desktop app needs to be open