# raDImod

Repository for _ab initio_ 3D modeling. Using information extraced from [ArchDBmap](https://github.com/jaumebonet/archdbmap) and [RADI](https://github.com/structuralbioinformatics/RADI) and passing it into Modeller this set of scripts allow to create protein models just form sequence information.

## Brief explanation of modules

This repository is built around three different modules.

### Builder

This module needs the alignment file obtained when executing **ArchDBmap**. _XXXX_X.archsearch_ file. Execution will output three different sets of files:
  - Alignment file to pass to `Models` module. _XXXX_X_ali.pir_
  - Alignment file to clean reference pdb for RMSD computations in the `Evaluation` module. _XXXX_X_dummy_ali.pir_
  - Pdb templates to pass to `Models` module.

 ### Models

 This module will create the three dimensional models using the alignment file _XXXX_X_ali.pir_ and the contacts file given by **RADI**. The latter must have _.out_ extension.

 ### Evaluation

 This module will assess the previously obtained models by ranking all the models by [DOPE score](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2242414/) along with the RMSD computed against the reference structure.


 ## Usage

 ### Requirements

  - argparse
  ```
  pip install argparse
  ```

  - Biopython
  ```
  pip install biopython
  ```

  - Modeller
    - If under a regular Python distribution see: https://salilab.org/modeller/9.24/release.html#deb
    - If using conda:
    ```
    conda config --add channels salilab
    conda install modeller
    ```
  ### Execution

  Basic execution of the module obtained with the `-h --help` option. Please, execute the scripts only with Python 2 (>= 2.7).

  #### Builder
  ```{console}
  usage: alignmentBuilder.py [-h] -i INPUT

  optional arguments:
    -h, --help            show this help message and exit
    -i INPUT, --input INPUT
                          path to folder with ArchDBmap output.

  ```

  #### Models
  ```{console}
  usage: structBuilder.py [-h] -i ARC -a ALI -s REALIGN -r RADI -p PDB
                        [-m MODELS]

  optional arguments:
    -h, --help            show this help message and exit
    -i ARC, --arc ARC     path to ArchDBmap output.
    -a ALI, --ali ALI     path to alignment file (ali.pir).
    -s REALIGN, --realign REALIGN
                          path to secondary structure prediction (.realign).
    -r RADI, --radi RADI  path to folder containing raDI output.
    -p PDB, --pdb PDB     path to folder containing template pdb structures.
    -m MODELS, --models MODELS
                          number of models to build.
  ```

  #### Evaluation
  ```{console}
  usage: modelEval.py [-h] -i INPUT -p PATH_ALI -c CODE

  optional arguments:
    -h, --help            show this help message and exit
    -i INPUT, --input INPUT
                          path to generated PDB structures.
    -p PATH_ALI, --path_ali PATH_ALI
                          path to dummy alignment.
    -c CODE, --code CODE  Code.

  ```
