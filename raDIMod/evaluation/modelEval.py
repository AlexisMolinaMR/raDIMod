import os
import argparse as ap
import glob
from modeller import *
from modeller.scripts import complete_pdb
from modeller import *
from modeller.automodel import *
from Bio.PDB import *
from Bio.PDB import PDBIO


def parseArg():

    parser = ap.ArgumentParser()

    parser.add_argument("-i", "--input", required=True, type=str, help="path to generated PDB structures.")
    parser.add_argument("-p", "--path_ali", required=True, type=str, help="path to dummy alignment.")
    parser.add_argument("-c", "--code", required=True, type=str, help="Code.")

    args = parser.parse_args()

    path = args.input
    path_ali = args.path_ali
    code = args.code

    return path, path_ali, code

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

def get_reference_pdb(path_ali, code):

    chain_A = {}
    chain_B = {}

    io = PDBIO()
    pdb_list = PDBList()


    pdb_list.retrieve_pdb_file(code, file_format="pdb", pdir=path_ali)

    with open(path_ali + "pdb" + code.lower() + ".ent", 'r') as pdb:
        for lines in pdb:
            line = lines.split()
            if line[0] == "ATOM" or (line[0] == "HETATM" and line[3] in ["MSE","KCX","LYS"]):
                if line[4] == "A":
                    if not chain_A:
                        chain_A["A"] = lines
                    else:
                        chain_A["A"] += lines

                elif line[4] == "B":
                    if not chain_B:
                        chain_B["B"] = lines
                    else:
                        chain_B["B"] += lines


    try:
        with open(path_ali + code + "_A.pdb", "w") as pdb_out_A:
            pdb_out_A.write(chain_A["A"])
        with open(path_ali + code + "_B.pdb", "w") as pdb_out_B:
            pdb_out_B.write(chain_B["B"])
    except:
        pass


def pdb_cleaner(best_struct, path_ali, code):

    dir = os.getcwd()
    os.chdir(path_ali)

    global chain_to_use
    chain_to_use = best_struct[1].split("/")[-1].split(".")[0].split("_")[1]

    log.verbose()
    env = environ()

    env.io.hetatom= True
    io_data.hetatm = True

    env.io.atom_files_directory = [path_ali]

    a = automodel(env,
                  alnfile = path_ali + best_struct[1].split("/")[-1].split(".")[0] + "_dummy_ali.pir",
                  knowns = code + "_" + chain_to_use,
                  sequence = best_struct[1].split("/")[-1].split(".")[0].split("_")[0])

    a.starting_model = 1
    a.ending_model = 1

    a.make()

    os.chdir(path_ali)


def rmsd(path, path_ali, ranked_dope_scores, code):

    global evaluation
    evaluation = []
    global ranked_evaluation
    ranked_evaluation = []

    reference_structure = path_ali +  code + ".B99990001.pdb"

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
    with open(path + best_struct[1].split("/")[-1].split(".")[0] + '_evaluation.txt', 'w') as dope_file:
        print("\nWriting strucutres ranking...\n")
        dope_file.write("Structure\t Dope score\t RMSD\n")
        for score_file in ranked_evaluation:
            if best_struct[1] == score_file[1]:
                dope_file.write(str(score_file[1].split('/')[-1]) + '\t ' + str(score_file[0]) + '\t ' + str(score_file[2]) + "\n")
            else:
                dope_file.write(str(score_file[1].split('/')[-1]) + '\t ' + str(score_file[0]) + '\t ' + str(score_file[2]) + '\n')

    print("DONE!")


def main():

    path, path_ali, code = parseArg()

    get_reference_pdb(path_ali=path_ali, code=code)
    dope_score(path=path)
    pdb_cleaner(best_struct=best_struct, path_ali=path_ali, code=code)
    rmsd(path=path, path_ali=path_ali, ranked_dope_scores=ranked_dope_scores, code=code)
    score_writer(path=path, ranked_evaluation=ranked_evaluation)


if __name__ == "__main__":
    main()
