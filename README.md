# RNA_LEMARCHAND_POSTIC

The 'training.py' script computes interactominc distances and pseudo-energy (scores) of known structures from pdb files.
To launch the script as a command line use the following command : python3 training.py input_file
The input file must be a pdb file, example with the test file : python3 training.py data/2jyh.pdb
The 'scoring.r' plots interaction profiles withthe score as a function of the distance.

The 'data' directory contains a pdb file used for testing and a file containing scores for each interactions.

The 'plot' directory contains all plots obtained with the scoring script.
