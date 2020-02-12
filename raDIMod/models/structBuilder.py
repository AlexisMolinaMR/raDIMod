# Load standard Modeller classes
from modeller import *
# Load the automodel class
from modeller.automodel import *
import argparse as ap
import glob
import os

def parseArg():

    parser = ap.ArgumentParser()

    parser.add_argument("-i", "--arc", required=True, type=str, help="path to ArchDBmap output.")
    parser.add_argument("-a", "--ali", required=True, type=str, help="path to alignment file (ali.pir).")
    parser.add_argument("-s", "--realign", required=True, type=str, help="path to secondary structure prediction (Realign).")
    parser.add_argument("-r", "--radi", required=True, type=str, help="path to folder containing raDI output.")
    parser.add_argument("-p", "--pdb", required=True, type=str, help="path to folder containing template pdb structures.")
    parser.add_argument("-m", "--models", required=False, type=int, help="number of models to build. Default set to 100 models")

    args = parser.parse_args()

    path_arc = args.arc
    path_ali = args.ali
    path_realign = args.realign
    path_radi = args.radi
    path_pdb = args.pdb
    num_models = args.models


    return path_arc, path_ali, path_realign, path_radi, path_pdb, num_models

codes_arc = []

seq = ""

temp_codes = []

lengths = []
rel_pos = []
pair = ()

global MI_dis
MI_dis = []
global DI_dis
DI_dis = []
global pos
pos = ()

def retrieve_codes(path_arc):
    file = open(path_arc).readlines()
    print("Opening " + path_arc + '\n')
    for lines in file:
        k = lines.split()
        if "#" not in k[0][0] and k[0] not in codes_arc:
            print("Retrieving code: " + k[0])
            codes_arc.append(k[0])

    print("{} codes retrieved!\n".format(len(codes_arc)))
    return codes_arc


def retrieve_template_codes(path_arc):
    file = open(path_arc).readlines()
    print("Retrieving template codes...")
    for lines in file:
        k = lines.split()
        try:
            temp_codes.append([k[4],k[1]+"*"])
        except:
            pass

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

    return tuple(temp_codes)

def secondary_structure_alpha(path_realign):

    s = ""
    fp = open(path_realign)
    for i, line in enumerate(fp):
        if i == 3:
            s = line
    #print(s)
    l = []
    for i in range(len(s)):
        if s[i] == "H" and s[i+1] == "H":
            l.append(i+1)
        elif s[i] == "H" and s[i+1] != "H":
            l.append(i+1)

    global posalpha
    posalpha = []
    posalpha.append(s.index("H")+1)

    i = 0
    while i < len(l)-1:
        if l[i] + 1 != l[i+1]:
            posalpha.append(l[i])
            posalpha.append(l[i+1])
        i += 1

    posalpha.append(len(s)-s[::-1].index("H"))

    #print(posalpha)
    return posalpha

def secondary_structure_beta(path_realign):
    s = ""
    fp = open(path_realign)
    for i, line in enumerate(fp):
        if i == 3:
            s = line
    l = []
    for i in range(len(s)):
    	if s[i] == "E" and s[i+1] == "E":
    		l.append(i+1)
    	elif s[i] == "E" and s[i+1] != "E":
    		l.append(i+1)

    global posbeta
    posbeta = []
    posbeta.append(s.index("E"))

    i = 0
    while i < len(l)-1:
    	if l[i] + 1 != l[i+1]:
    		posbeta.append(l[i])
    		posbeta.append(l[i+1])
    	i += 1

    posbeta.append(len(s)-s[::-1].index("E"))

    return posbeta
    #print(posbeta)

def raDI_contacts(path_radi):
    for filename in glob.glob(os.path.join(path_radi, '*.out')):
        file = open(filename).readlines()
        for lines in file:
            k = lines.split()
            try:
                if (k[19] == "yes" or k[19] == "no") and (k[1] != "-" or k[2] != "-"):
                    pos = (k[0],k[1],k[2],k[7],k[8],'MI')
                    MI_dis.append(pos)
                    pos = ()

                if (k[20] == "yes" or k[20] == "no") and (k[1] != "-" or k[2] != "-"):
                    pos = (k[0],k[4],k[5],k[9],k[10],'DI')
                    DI_dis.append(pos)
                    pos = ()
            except:
                pass

    return MI_dis, DI_dis

class MyModel(loopmodel):
    def special_restraints(self,alnfile):
        rsr = self.restraints
        at = self.atoms

        #addition of distance restraints. THIS SHOULD BE DONE 5 TIMES, one for each alphabet.
        for c in range(0,160):
            rsr.add(forms.gaussian(group=physical.xy_distance,feature=features.distance(at['CA: {}'.format(MI_dis[c][1])],at['CA: {}'.format(MI_dis[c][2])]), mean = float(MI_dis[c][3]), stdev = 2.5))

        for d in range(0,160):
            rsr.add(forms.gaussian(group=physical.xy_distance,feature=features.distance(at['CA: {}'.format(DI_dis[d][1])],at['CA: {}'.format(DI_dis[d][2])]), mean = float(DI_dis[d][3]), stdev=2.5))

        for w in range(0,len(posalpha),2):
            rsr.add(secondary_structure.alpha(self.residue_range('{}:'.format(posalpha[w]), '{}:'.format(posalpha[w+1]))))
        for r in range(0,len(posbeta),2):
            rsr.add(secondary_structure.strand(self.residue_range('{}:'.format(posbeta[r]), '{}:'.format(posbeta[r+1]))))


def main():
    path_arc, path_ali, path_realign, path_radi, path_pdb, num_models = parseArg()

    retrieve_codes(path_arc=path_arc)
    retrieve_template_codes(path_arc=path_arc)
    retrieve_relative_postions(temp_codes=temp_codes)
    fill_temp_codes_tuple(temp_codes=temp_codes, rel_pos=rel_pos)

    secondary_structure_alpha(path_realign=path_realign)
    secondary_structure_beta(path_realign=path_realign)

    raDI_contacts(path_radi=path_radi)

    log.verbose() # request verbose output
    env = environ() # create a new MODELLER environment to build this model in

    # declare directories for input atom files
    env.io.atom_files_directory = ['.', path_pdb]

    a = MyModel(env,
                alnfile = path_ali,
                knowns = ('4KD5_C_102', '4KD5_C_115', '4KD5_C_119', '4KD5_C_150', '4KD5_C_181', '4KD5_C_187', '4KD5_C_201', '4KD5_C_51', '4KD5_C_71', '4KD5_C_80', '4KD5_C_95', '2H5Y_B_11', '2H5Y_B_115', '2H5Y_B_124', '2H5Y_B_14', '2H5Y_B_146', '2H5Y_B_171', '2H5Y_B_219', '2H5Y_B_32', '2H5Y_B_39', '2H5Y_B_5', '2H5Y_B_54', '2H5Y_B_60', '2H5Y_B_74'),
                sequence = codes_arc[0],
                assess_methods = assess.DOPE )

    a.md_level=refine.slow
    a.loop.starting_model=1
    a.loop.ending_model=4
    a.loop.md_level=refine.slow

    a.starting_model= 0             # index of the first model
    a.ending_model  = num_models    # index of the last model
                                    # (determines how many models to calculate)
    a.make()                        # do the actual homology modeling


if __name__ == "__main__":
    main()
