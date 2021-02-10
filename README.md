# X-EISD: Extended Experimental Inferential Structure Determination
X-EISD is a Bayeian approach to perform Experimental Inferential Structure Determination of ensembles for intrinsically disordered proteins.

Program version: 0.2.0  -  Latest update: April 03, 2020


## Installation:
You can install eisd from the source:

    git clone https://github.com/THGLab/X-EISD.git
    cd X-EISD
    pip install -e .

## Git LFS data quota limit
If you are not able to download the data using git lfs, please try to download it directly from the following link: 
https://datadryad.org/stash/share/yrRQFe-bjpDqupDFtT_pxYJjrbT6cdhhUGNTVl-JLSA

Please unzip and copy the contents to the `data` directory in this repository. The size of the compressed data is 1.14 GB. 

## Dependencies:
The dependencies for X-EISD are only numpy and pandas libraries. A proper installation of the package will also install
the requirements.

## Citation:
Please cite the use of X-EISD as:


    (1) James Lincoff, Mojtaba Haghighatlari, Mickael Krzeminski, Joao Teixeira, Gregory Neal-Gomes, Claudu Gardinaru, Julie Forman-Kay, Teresa Head-Gordon, https://www.nature.com/articles/s42004-020-0323-0
    (2) David H. Brookes, and Teresa Head-Gordon, Experimental Inferential Structure Determination of Ensembles for Intrinsically Disordered Proteins, JACS 138, 2016, 4530-4538 DOI: 10.1021/jacs.6b00351

## Getting Started 
You can either follow along the sample_script.py in the repository or use the commnd line interface to run eisd: 

    python sample_script.py     # first modify sample_script based upon your request 
      
or

    eisdshell -d/--datapath, -m/--mode, -s/--structure, -e/--epochs, -o/--output



