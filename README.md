# X-EISD: Extended Experimental Inferential Structure Determination
X-EISD is a Bayeian approach to perform Experimental Inferential Structure Determination of ensembles for intrinsically disordered proteins.

Program version: 0.2.0  -  Latest update: April 03, 2020


## Installation:
You can install eisd from the source:

    git clone https://github.com/THGLab/X-EISD.git
    cd X-EISD
    pip install -e .

## Dependencies:
The dependencies for X-EISD are only numpy and pandas libraries. A proper installation of the package will also install
the requirements.

## Citation:
Please cite the use of X-EISD as:


    (1) James Lincoff, Mickael Krzeminski, Mojtaba Haghighatlari, Joao Teixeira, Gregory Neal-Gomes, Claudu Gardinaru, Julie Forman-Kay, Teresa Head-Gordon, arXiv:1912.12582, 2019
    (2) David H. Brookes, and Teresa Head-Gordon, Experimental Inferential Structure Determination of Ensembles for Intrinsically Disordered Proteins, JACS 138, 2016, 4530-4538 DOI: 10.1021/jacs.6b00351

## Getting Started 
You can either follow along the sample_script.py in the repository or use the commnd line interface to run eisd: 

    python sample_script.py     # first modify sample_script based upon your request 
      
or

    eisdshell -d/--datapath, -m/--mode, -s/--structure, -e/--epochs, -o/--output



