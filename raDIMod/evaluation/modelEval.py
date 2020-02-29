import os
import argparse as ap
import glob
from modeller import *
from modeller.scripts import complete_pdb
from modeller import *
from modeller.automodel import *
from Bio.PDB import *


def parseArg():

    parser = ap.ArgumentParser()

    parser.add_argument("-i", "--input", required=True, type=str, help="path to generated PDB structures.")
    parser.add_argument("-r", "--reference", required=True, type=str, help="path to reference PDB structure.")
    parser.add_argument("-p", "--path_ali", required=True, type=str, help="path to dummy alignment.")

    args = parser.parse_args()

    path = args.input
    reference_path = args.reference
    path_ali = args.path_ali

    return path, reference_path, path_ali

def dope_score(path):

    global dope_scores
    dope_scores = []
    global ranked_dope_scores
    ranked_dope_scores = []
    global best_struct
    best_struct = [-0.1,"dummy"]

    env = environ()
    env.libs.topology.read(file='$(LIB)/top_heav.lib')
    env.libs.parameters.read(file='$(LIB)/par.lib')

    for filename in glob.glob(os.path.join(path, '*.pdb')):

        print("Opening " + filename)

        mdl = complete_pdb(env, filename)

        atmsel = selection(mdl)

        print("Computing DOPE score...\n")
        score = atmsel.assess_dope()

        dope_scores.append((score,filename))

        if best_struct[0] > score:

            best_struct[0] = score
            best_struct[1] = filename

    ranked_dope_scores = sorted(dope_scores, key=lambda x: x[0])

    return ranked_dope_scores, best_struct

def pdb_cleaner(best_struct, reference_path, path_ali):

    log.verbose()
    env = environ()

    env.io.atom_files_directory = [reference_path] #get only path

    a = automodel(env,
                  alnfile = path_ali + best_struct[1].split("/")[-1].split(".")[0] + "_dummy_ali.pir",
                  knowns = ("1ATG"),
                  sequence = best_struct[1].split("/")[-1].split(".")[0])

    a.starting_model = 1
    a.ending_model = 1

    a.make()

#    ref = a.get_model_filename()

#    return ref

def rmsd(path, ranked_dope_scores):

    global evaluation
    evaluation = []
    global ranked_evaluation
    ranked_evaluation = []

    reference_structure = "/home/alexis/Desktop/raDIMod/raDIMod/evaluation/1ATG_A.B99990001.pdb"

    parser = PDBParser()

    reference_struct = parser.get_structure('reference', reference_structure)
    reference_struct = list(reference_struct.get_atoms())

    for filename in glob.glob(os.path.join(path, '*.pdb')):
        structure = parser.get_structure('struct', filename)
        structure = list(structure.get_atoms())

        print("Superimposing structures...\n")

        sup = Superimposer()

        print("Computing RMSD...\n")

        sup.set_atoms(reference_struct, structure)

        for score_file in ranked_dope_scores:
            if filename.split("/")[-1]  == score_file[1].split('/')[-1]:
                score_file = score_file + (sup.rms,)
                evaluation.append(score_file)

        print("RMSD between reference structure and " + filename.split("/")[-1] + " is {}".format(sup.rms))

    ranked_evaluation = sorted(evaluation, key=lambda x: x[0])

    return ranked_evaluation

def score_writer(path, ranked_evaluation):
    with open(path + 'evaluation.txt', 'w') as dope_file:
        print("\nWriting strucutres ranking...\n")
        dope_file.write("Structure\t Dope score\t RMSD\n")
        for score_file in ranked_evaluation:
            if best_struct[1] == score_file[1]:
                dope_file.write(str(score_file[1].split('/')[-1]) + '\t ' + str(score_file[0]) + '\t ' + str(score_file[2]) + "\n")
            else:
                dope_file.write(str(score_file[1].split('/')[-1]) + '\t ' + str(score_file[0]) + '\t ' + str(score_file[2]) + '\n')

    print("DONE!")


def main():

    path, reference_path, path_ali = parseArg()

    dope_score(path=path)
    pdb_cleaner(best_struct=best_struct, reference_path=reference_path, path_ali=path_ali)
    rmsd(path=path, ranked_dope_scores=ranked_dope_scores)
    score_writer(path=path, ranked_evaluation=ranked_evaluation)


if __name__ == "__main__":
    main()
