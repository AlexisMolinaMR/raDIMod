import argparse as ap
import os
import sys
import glob
from Bio.PDB import *
from Bio.PDB import PDBIO

def parseArg():

    parser = ap.ArgumentParser()

    parser.add_argument("-i", "--input", required=True, type=str, help="path to ArchDBmap output.")

    args = parser.parse_args()

    path = args.input

    return path


codes_arc = []

seq = ""

temp_codes = []

lengths = []
rel_pos = []
pair = ()

def retrieve_codes(path):

    txt_files = [f for f in os.listdir(path) if f.endswith('.archsearch')]
    if len(txt_files) != 1:
        raise ValueError('should be only one txt file in the current directory')

    global filename
    filename = txt_files[0]
    file = open(path + filename).readlines()
    print("Opening " + path + '\n')

    for lines in file:
        k = lines.split()
        if len(k) > 0:
            if "#" not in k[0][0] and k[0] not in codes_arc:
                print("Retrieving code: " + k[0])
                codes_arc.append(k[0])


    print("{} codes retrieved!\n".format(len(codes_arc)))
    return codes_arc, filename

def retrieve_sequence(path,codes_arc,filename):
    file = open(path + filename).readlines()
    print("Retrieving sequence...")
    for lines in file:
        k = lines.split()
        if len(k) > 0:
            if k[0] == codes_arc[0]:
                try:
                    seq = k[1]
                except:
                    pass

    if len(seq) > 0:
        print("Sequence retrieved successfully!\n")
        return seq
    else:
        raise Exception("Couldn't find sequence...take a look at ArchDBmap output.")

def retrieve_template_codes(path, codes_arc, filename):
    file = open(path + filename).readlines()
    print("Retrieving template codes...")
    for lines in file:
        k = lines.split()
        if len(k) > 0 and k[0][0] != '#' and k[0] != codes_arc[0]:
            try:
                temp_codes.append([k[4],k[1]+"*"])
            except:
                pass

        elif len(k) == 0:
            if len(temp_codes) > 0:
                print("Codes retrieved successfully!\n")
                return temp_codes
            else:
                raise Exception("Couldn't find template codes...take a look at ArchDBmap output.")


def retrieve_relative_postions(temp_codes):
    print("Retrieving relative postions...")
    for i in temp_codes:
        lengths = i[1].split(".")
        for j in (lengths):
            if len(j) > 1:
                pair = (i[0].split("_")[2],str(int(i[0].split("_")[2])+len(j)),i[0].split("_")[1])
                rel_pos.append(pair)

    if len(rel_pos) > 0:
        print("Relative positons retrieved successfully!\n")
        return rel_pos
    else:
        raise Exception("Couldn't extract relative positions...take a look at ArchDBmap output.")


def fill_temp_codes_tuple(temp_codes, rel_pos):
    print("Pairing template codes and positions...")
    for ele in range(len(rel_pos)):
        temp_codes[ele].append(rel_pos[ele])

    for i in range(len(temp_codes)):
        temp_codes[i][1] = temp_codes[i][1].replace(".","-")

    print("Done!\n")

    return temp_codes


def aliBuild(codes_arc, seq, temp_codes, path):

    print("Building alignment...")
    print(temp_codes)
    seq += "*"
    with open(path + codes_arc[0] + "_ali.pir", "w") as f1:

        f1.write(">P1;{}\n".format(codes_arc[0]))
        f1.write("SEQUENCE:{}: : : : : : : SEQUENCE :(SEQU) .\n".format(codes_arc[0]))
        f1.write(seq+'\n')
        for i in temp_codes:
            f1.write(">P1;{}\n".format(i[0]))
            f1.write("structureX:{}:{}:{}:{}:::::\n".format(i[0],i[2][0],i[2][2],i[2][1]))
            f1.write(i[1]+'\n')


    print(codes_arc[0] + "_ali.pir" + " built successfully!\n")

def dummy_alibuilder(path, seq, codes_arc):

    first_residues = 0

    pdb_list = PDBList()
    pdb_list.retrieve_pdb_file(codes_arc[0].split("_")[0], file_format="pdb", pdir=path)

    with open(path + "pdb" + codes_arc[0].split("_")[0].lower() + ".ent", 'r') as pdb:
        for lines in pdb:
            line = lines.split()
            if line[0] == "ATOM" or (line[0] == "HETATM" and line[3] in ["MSE","KCX","LYS"]):
                if line[4] == codes_arc[0].split("_")[1]:
                    first_residues = line[5]
                    break


    seq += "*"

    with open(path + codes_arc[0] + "_dummy_ali.pir", 'w') as dum_ali:
        print("Creating dummy alignment for model evaluation...")
        dum_ali.write(">P1;{}\n".format(codes_arc[0].split("_")[0]))
        dum_ali.write("SEQUENCE:{}: : : : : : : SEQUENCE :(SEQU) .\n".format(codes_arc[0].split("_")[0]))
        dum_ali.write(seq+'\n')
        dum_ali.write(">P1;{}\n".format(codes_arc[0]))
        dum_ali.write("structureX:{}:{}:{}:{}:::::\n".format(codes_arc[0],first_residues,codes_arc[0].split("_")[1],len(seq[:-1])))
        dum_ali.write(seq+'\n')

    print("Dummy alignment created!\n")

def pdb_cutter(path, temp_codes):

    amino = {'C': 'CYS', 'D': 'ASP', 'S': 'SER', 'Q': 'GLN', 'K': 'LYS',
             'I': 'ILE', 'P': 'PRO', 'T': 'THR', 'F': 'PHE', 'N': 'ASN',
             'G': 'GLY', 'H': 'HIS', 'L': 'LEU', 'R': 'ARG', 'W': 'TRP',
             'A': 'ALA', 'V':'VAL', 'E': 'GLU', 'Y': 'TYR', 'M': 'MET'}

    pdb_list = PDBList()

    temp_store = []
    for i in temp_codes:
        pdb_list.retrieve_pdb_file(i[0].split("_")[0], file_format="pdb", pdir=path)


    for i in temp_codes:
        count_pos = 0
        frag_store = []
        all_letter = False

        for pos in range(len(i[1])):


            if i[1][pos] not in "-*" and len(temp_store) == 0 and len(frag_store) == 0:  #first fragment first res
                frag_store.append(amino[i[1][pos]])
                frag_store.append(i[0].split("_")[2])
                frag_store.append(i[0].split("_")[1])

            elif i[1][pos] not in "-*" and len(temp_store) == 0 and len(frag_store) != 0: # first fragment
                count_pos += 1
                frag_store.append(amino[i[1][pos]])
                frag_store.append(int(i[0].split("_")[2]) + count_pos)
                frag_store.append(i[0].split("_")[1])

            elif i[1][pos] not in "-*" and len(temp_store) != 0:
                frag_store.append(amino[i[1][pos]])
                frag_store.append(int(i[0].split("_")[2]) + count_pos)
                frag_store.append(i[0].split("_")[1])
                count_pos += 1

            elif i[1][pos] == "-" and len(frag_store) != 0 and not all_letter:
                rest_seq = ""
                all_seq = ""

                with open(path + "pdb" + i[0].split("_")[0].lower() + ".ent", 'r') as pdb:
                    for lines in pdb:
                        line = lines.split()
                        if line[0] == "ATOM" and line[1] == '1':
                            first_residu_pos = int(line[5])
                            break

                for ele in range(len(i[1][pos:])):
                    if i[1][pos:][ele] not in "-*":
                        rest_seq += i[1][pos:][ele]
                        all_letter = True


                if all_letter:
                    p = PDBParser()
                    structure = p.get_structure(i[0].split("_")[0].lower(), path + "pdb" + i[0].split("_")[0].lower() + ".ent")
                    ppb = CaPPBuilder()

                    polypeptide = ppb.build_peptides(structure)[0]
                    all_seq = polypeptide.get_sequence()


                if len(rest_seq) != 0:
                    for every_res in range(0, len(all_seq)):
                        if rest_seq == all_seq[every_res:every_res+len(rest_seq)]:
                            for posi in range(len(rest_seq)):
                                frag_store.append(amino[rest_seq[posi]])
                                frag_store.append(first_residu_pos + every_res + posi)
                                frag_store.append(i[0].split("_")[1])

                break

        frag_store.append(i[0])
        temp_store.append(frag_store)


    for frag in temp_store:
        f = open(path + frag[-1] + ".pdb", 'w')
        f.close()
        for res in range(0, len(frag)-1, 3):
            with open(path + "pdb" + frag[-1].split("_")[0].lower() + ".ent", 'r') as pdb:
                for lines in pdb:
                    line = lines.split()
                    if line[0] == "ATOM" and line[3] == frag[res] and line[4] == str(frag[res+2]) and line[5] == str(frag[res+1]):
                        with open(path + frag[-1] + ".pdb", 'a') as pdb_out:
                            pdb_out.write(lines)



    print("Templates writen and stored in " + path)


def main():

    path = parseArg()

    retrieve_codes(path=path)
    retrieve_sequence(path=path,codes_arc=codes_arc, filename=retrieve_codes(path=path)[1])
    retrieve_template_codes(path=path, codes_arc=codes_arc, filename=retrieve_codes(path=path)[1])
    retrieve_relative_postions(temp_codes=temp_codes)
    fill_temp_codes_tuple(temp_codes=temp_codes, rel_pos=rel_pos)
    aliBuild(codes_arc=codes_arc,seq=retrieve_sequence(path=path,codes_arc=codes_arc,filename=retrieve_codes(path=path)[1]),temp_codes=temp_codes, path=path)
    dummy_alibuilder(path=path, seq=retrieve_sequence(path=path,codes_arc=codes_arc,filename=retrieve_codes(path=path)[1]), codes_arc=codes_arc)
    pdb_cutter(path=path, temp_codes=temp_codes)

    print("\nAhead with your models!")


if __name__ == "__main__":
    main()
