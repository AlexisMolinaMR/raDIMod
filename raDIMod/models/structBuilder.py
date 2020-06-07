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
    parser.add_argument("-s", "--realign", required=True, type=str, help="path to secondary structure prediction (.realign).")
    parser.add_argument("-r", "--radi", required=True, type=str, help="path to folder containing raDI output.")
    parser.add_argument("-p", "--pdb", required=True, type=str, help="path to folder containing template pdb structures.")
    parser.add_argument("-m", "--models", required=False, type=int, help="number of models to build.")

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
        try:
            if "#" not in k[0][0] and k[0] not in codes_arc:
                print("Retrieving code: " + k[0])
                codes_arc.append(k[0])
        except:
            pass
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
                try:
                    pair = (i[0].split("_")[2],str(int(i[0].split("_")[2])+len(j)),i[0].split("_")[1])
                    rel_pos.append(pair)
                except:
                    pass

    if len(rel_pos) > 0:
        print("Relative positons retrieved successfully!\n")
        return rel_pos
    else:
        raise Exception("Couldn't extract relative positions...take a look at ArchDBmap output.")


def fill_temp_codes_tuple(temp_codes, rel_pos):
    print("Pairing template codes and positions...\n")
    for ele in range(len(rel_pos)):
        temp_codes[ele].append(rel_pos[ele])

    for i in range(len(temp_codes)):
        temp_codes[i][1] = temp_codes[i][1].replace(".","-")

    print("Done!\n")

    global cod
    cod = ()

    for tem in temp_codes:
        print(tem[0])
        cod += (tem[0],)

    cod = cod[1:]

    return cod

def secondary_structure_alpha(path_realign):

    global there_is_alpha
    there_is_alpha = True

    s = ""
    print("Opening " + path_realign + '\n')
    fp = open(path_realign)
    for i, line in enumerate(fp):
        if i == 3:
            s = line

    l = []

    if "H" in s:
        print("Retrieving alpha helix positions...")
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


        if len(posalpha) > 0:
            print("Alpha positions retrieved!")
            return posalpha
        else:
            raise Exception("Unable to retrieve alpha positions. Take a look at the realign file.")
    else:
        "No alpha positions in SS prediction. (.realign)"
        there_is_alpha = False

    return there_is_alpha

def secondary_structure_beta(path_realign):

    global there_is_beta
    there_is_beta = True

    s = ""
    fp = open(path_realign)
    for i, line in enumerate(fp):
        if i == 3:
            s = line
    l = []

    if "E" in s:
        print("Retrieving beta sheet positions...\n")

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

        if len(posbeta) > 0:
            print("Beta positions retrieved!\n")
            return posbeta
        else:
            raise Exception("Unable to retrieve beta positions. Take a look at the realign file.")
    else:
        "No beta positions in SS prediction. (.realign)"
        there_is_beta = False

    return there_is_beta

def raDI_contacts(path_radi):
    for filename in glob.glob(os.path.join(path_radi, '*.out')):
        file = open(filename).readlines()
        print("Opening" + filename + '\n')
        for lines in file:
            k = lines.split()
            try:
                if (k[19] == "yes" or k[19] == "no") and (k[1] != "-" or k[2] != "-"):
                    pos = (k[0],k[1],k[2],k[7],k[8],'MI')
                    MI_dis.append(pos)
                    print("Retrieving Mutual Information...")
                    pos = ()

                if (k[20] == "yes" or k[20] == "no") and (k[1] != "-" or k[2] != "-"):
                    pos = (k[0],k[4],k[5],k[9],k[10],'DI')
                    DI_dis.append(pos)
                    print("Retrieving Direct Information")
                    pos = ()
            except:
                pass

    if len(MI_dis) > 0 and len(DI_dis) > 0:
        print("Mutual and Direct information retrieved!\n")
        return MI_dis.sort(key = lambda x: int(x[0])), DI_dis.sort(key = lambda x: int(x[0]))
    else:
        raise Exception("Couldn't retrieve MI and/or DI...take a look at raDI output.")

class MyModel(loopmodel):
    print("Calling MODELLER...")
    def special_restraints(self,alnfile):
        print("Adding special restraints...")
        rsr = self.restraints
        at = self.atoms

        for c in range(0,160):
            rsr.add(forms.gaussian(group=physical.xy_distance,
                                    feature=features.distance(at['CA: {}'.format(MI_dis[c][1])],
                                                                at['CA: {}'.format(MI_dis[c][2])]),
                                                                    mean = float(MI_dis[c][3]), stdev = 2.5))

        for d in range(0,160):
            rsr.add(forms.gaussian(group=physical.xy_distance,
                                    feature=features.distance(at['CA: {}'.format(DI_dis[d][1])],
                                                                at['CA: {}'.format(DI_dis[d][2])]),
                                                                    mean = float(DI_dis[d][3]), stdev=2.5))
        if there_is_alpha:
            for w in range(0,len(posalpha),2):
                rsr.add(secondary_structure.alpha(self.residue_range('{}:'.format(posalpha[w]), '{}:'.format(posalpha[w+1]))))
        if there_is_beta:
            for r in range(0,len(posbeta),2):
                rsr.add(secondary_structure.strand(self.residue_range('{}:'.format(posbeta[r]), '{}:'.format(posbeta[r+1]))))

        print("To verbose mode ->\n")

def main():
    path_arc, path_ali, path_realign, path_radi, path_pdb, num_models = parseArg()

    retrieve_codes(path_arc=path_arc)
    retrieve_template_codes(path_arc=path_arc)
    retrieve_relative_postions(temp_codes=temp_codes)
    fill_temp_codes_tuple(temp_codes=temp_codes, rel_pos=rel_pos)

    secondary_structure_alpha(path_realign=path_realign)
    secondary_structure_beta(path_realign=path_realign)

    raDI_contacts(path_radi=path_radi)

    log.verbose()
    env = environ()

    env.io.atom_files_directory = ['.', path_pdb]

    a = MyModel(env,
                alnfile = path_ali,
                knowns = cod,
                sequence = codes_arc[0],
                assess_methods = assess.DOPE )

    a.md_level=refine.slow
    a.loop.starting_model=1
    a.loop.ending_model=4
    a.loop.md_level=refine.slow

    a.starting_model= 0
    a.ending_model  = num_models

    a.make()


if __name__ == "__main__":
    main()
