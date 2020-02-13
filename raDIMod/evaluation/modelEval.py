import os
from modeller import *
from modeller.scripts import complete_pdb
import argparse as ap
import glob


'''
def rmsd(best_struct, reference):

    reference_structure = os.path.basename(reference)
    best_model = os.path.basename(best_struct[0][1])

    parser = PDBParser()
    structure = parser.get_structure('best', best_struct[1])
    reference_struct = parser.get_structure('reference', reference)

    reference_struct = list(reference_struct.get_atoms())
    structure = list(structure.get_atoms())

    sup = Superimposer()
    sup.set_atoms(reference_struct, structure)

    print(sup.rms)

    reader = Reader().readThisFile("my_trajectory.pdb").gettingOnlyCAs()
'''


def parseArg():

    parser = ap.ArgumentParser()

    parser.add_argument("-i", "--input", required=True, type=str, help="path to generated PDB structures.")
#    parser.add_argument("-r", "--reference", required=True, type=str, help="path to reference PDB structure.")

    args = parser.parse_args()

    path = args.input
#    reference = args.reference

    return path

def dope_score(path):

    global dope_scores
    dope_scores = []
    global best_struct
    best_struct = [-0.1,"dummy"]

    env = environ()
    env.libs.topology.read(file='$(LIB)/top_heav.lib')
    env.libs.parameters.read(file='$(LIB)/par.lib')

    for filename in glob.glob(os.path.join(path, '*.pdb')):
        mdl = complete_pdb(env, filename)

        atmsel = selection(mdl)

        score = atmsel.assess_dope()

        dope_scores.append((score,filename))

        if best_struct[0] > score:

            best_struct[0] = score
            best_struct[1] = filename


    return dope_scores, best_struct

def score_writer(path, dope_scores, best_struct):
    with open(path + 'dope_scores.txt', 'w') as dope_file:
        dope_file.write("Structure\t Dope score\n")
        for score_file in dope_scores:
            if best_struct[1] == score_file[1]:
                dope_file.write(str(score_file[1].split('/')[-1]) + '\t ' + str(score_file[0]) + "\tBest structure\n")
            else:
                dope_file.write(str(score_file[1].split('/')[-1]) + '\t ' + str(score_file[0]) + '\n')


def main():

    path = parseArg()

    dope_score(path=path)
    score_writer(path=path, dope_scores=dope_scores, best_struct=best_struct)

#    rmsd(best_struct=best_struct, reference=reference)

if __name__ == "__main__":
    main()
