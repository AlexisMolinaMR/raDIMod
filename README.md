# raDIMod

 Repository for building 3D models using [ArchDBmap](https://github.com/jaumebonet/archdbmap), [RADI](https://github.com/user/repo/blob/branch/other_file.md) and Modeller.

 ## Brief explanation of modules

This repository is build around three different independent modules.

### Builder

This module needs the alignment file obtained when executing **ArchDBmap**. _XXXX.archsearch_ file. 
Execution will output three different sets of things:
- Alignment file to pass to `Models` module. _XXXX_X_ali.pir_
- Alignment to clean reference pdb for RMSD computations in the `Evaluation` module.  _XXXX_X_dummy_ali.pir_
- Pdb templates to pass to `Models` module.

### Models

This module will create the models using the alignment file _XXXX_X_ali.pir_ and the contacts file given by **RADI** that must have _.out_ extension.

### Evaluation

This module will assess the previously obtained models by ranking all the models by [DOPE score](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2242414/) along with the RMSD with the reference structure.

## Usage
 
### Requirements
 
- argparse
 ```console
pip install argparse
 ```
- Biopython
 ```console
pip install biopython
 ```
- Modeller
  - If under a regular Python distribution see: https://salilab.org/modeller/9.24/release.html#deb
  - If using conda:
 ```console
conda config --add channels salilab
conda install modeller 
```
   
### Execution

Basic execution, for further explanation execute the module with the `-h --help` option.
 
#### Builder

 ```console
 python2 -i /home/user/XXXX/ 
 ```
 
 #### Models
 
  ```console
 python2 -i /home/user/XXXX/ 
 ```
 
 ##### Evaluation

 ```console
 python2 -i /home/user/XXXX/ 
 ```
 
