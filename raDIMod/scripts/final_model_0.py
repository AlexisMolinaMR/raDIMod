# Homology model by the automodel class

# Load standard Modeller classes
from modeller import *
# Load the automodel class
from modeller.automodel import *
# The module sys provides access to some variables used or maintained
# by the interpreter and to functions that interact strongly with the
# interpreter. Needed for creating the MODELLER environment below.
import sys
import argparse
import os

def parse_user_arguments(*args, **kwds):
    parser = argparse.ArgumentParser(description='Files needed for building the homology model')
    parser.add_argument('-p', '--Alignment file', dest = 'arcfile', action = 'store',
                        help = 'Alignment file for secodary structures (CODE.archsearch_sequence). It will be transformed to a .pir file.')
    parser.add_argument('-f', '--Realign file', dest = 'realign', action = 'store',
                        help = 'Alignment file in .fa format containing the secodary structure sequence.')
    options = parser.parse_args()
    return options
#    print(options)

#print(parse_user_arguments)

################################################################################

def main():

    options = parse_user_arguments()
#    print(options)
    arcfile = options.arcfile
#    print (alnfile)
    realign = options.realign
#    print (tuple(knowns))


    ### Piece of script for obtaining the alignment.pir file ###

    f = arcfile
    codes_arc = []
    seq = ""
    temp_codes = []

    file = open(f).readlines()
    for lines in file:
        k = lines.split()
        if "#" not in k[0][0] and k[0] not in codes_arc:
            codes_arc.append(k[0])
    file = open(f).readlines()
    for lines in file:
        k = lines.split()
        if k[0] == codes_arc[0]:
            try:
#                print(k[1])
                seq = k[1]
            except:
                pass
    file = open(f).readlines()
    for lines in file:
        k = lines.split()
        try:
            temp_codes.append([k[4],k[1]+"*"])
        except:
            pass

    #print(seq)
    lengths = []
    rel_pos = []
    pair = ()
    #print(temp_codes)
    for i in temp_codes:
        lengths = i[1].split(".")
        for j in (lengths):
            if len(j) > 1:
                pair = (i[0].split("_")[2],str(int(i[0].split("_")[2])+len(j)),i[0].split("_")[1])
                rel_pos.append(pair)


    #print(rel_pos)

    for ele in range(len(rel_pos)):
        temp_codes[ele].append(rel_pos[ele])

    for i in range(len(temp_codes)):
        temp_codes[i][1] = temp_codes[i][1].replace(".","-")

    seq += "*"
    f1 =  open("ali.pir", "w")

    f1.write(">P1;{}\n".format(codes_arc[0]))
    f1.write("SEQUENCE:{}: : : : : : : SEQUENCE :(SEQU) .\n".format(codes_arc[0]))
    f1.write(seq+'\n')
    for i in temp_codes:
        f1.write(">P1;{}\n".format(i[0]))
        f1.write("structureX:{}:{}:{}:{}:::::\n".format(i[0],i[2][0],i[2][2],i[2][1]))
        f1.write(i[1]+'\n')

    alnfile = "ali.pir"
    ################################################################################

    ################################################################################

    ### Piece of script for obtaining the positions of every alpha and beta region

    s = ""
    fp = open(realign)
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

    l = []
    for i in range(len(s)):
    	if s[i] == "E" and s[i+1] == "E":
    		l.append(i+1)
    	elif s[i] == "E" and s[i+1] != "E":
    		l.append(i+1)

    posbeta = []
    posbeta.append(s.index("E"))

    i = 0
    while i < len(l)-1:
    	if l[i] + 1 != l[i+1]:
    		posbeta.append(l[i])
    		posbeta.append(l[i+1])
    	i += 1

    posbeta.append(len(s)-s[::-1].index("E"))

    #print(posbeta)


    ################################################################################

    log.verbose() # request verbose output
    env = environ() # create a new MODELLER environment to build this model in

    # declare directories for input atom files
    env.io.atom_files_directory = ['.', '../atom_files']

    #retrieval of those first  MI (MI_dis) and DI (DI_dis) position we are sure about.
    # STILL TO BE DETERMINED BY FILE DISPOSITON

    MI_dis = []
    DI_dis = []
    pos = ()

    #program iterates through the 4 files containing the distances and retrives them by
    #appending them in a list, inside of a tuple.
    #FIRST element in the tuple is the ranked position.
    #SECOND and THIRD element are the postions of the pairs.
    #FOURTH element accounts for the distance between CA.
    #FIFTH element is the nearest MI or DI distance
    #The last entry is for distinguish between MUTUAL and DIRECT information.


    for z in range(0,4):
        f = "{}_raDI_filter9_ra{}_DI.out".format(codes_arc[0],z)
        file = open(f).readlines()
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

    MI_dis = sorted(MI_dis, key = lambda x: int(x[0]))

    DI_dis = sorted(DI_dis, key = lambda x: int(x[0]))

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

    knowns = []

    for codes in temp_codes:
        knowns.append(codes[0])

    print("TEMP_CODES:")
    print(temp_codes)


    knowns = tuple(knowns)
    print()
    print(knowns)

    a = MyModel(env,
                alnfile = "ali_mod.pir",
                knowns = ('4KD5_C_102', '4KD5_C_115', '4KD5_C_119', '4KD5_C_150', '4KD5_C_181', '4KD5_C_187', '4KD5_C_201', '4KD5_C_51', '4KD5_C_71', '4KD5_C_80', '4KD5_C_95', '2H5Y_B_11', '2H5Y_B_115', '2H5Y_B_124', '2H5Y_B_14', '2H5Y_B_146', '2H5Y_B_171', '2H5Y_B_219', '2H5Y_B_32', '2H5Y_B_39', '2H5Y_B_5', '2H5Y_B_54', '2H5Y_B_60', '2H5Y_B_74'),
                sequence = codes_arc[0],
                assess_methods = assess.DOPE )


    #Loopmodelling
    a.md_level=refine.slow
    a.loop.starting_model=1
    a.loop.ending_model=4
    a.loop.md_level=refine.slow

    a.starting_model= 0             # index of the first model
    a.ending_model  = 0	            # index of the last model
                                    # (determines how many models to calculate)
    a.make()                        # do the actual homology modeling

    #DOPE energies ranking
    # Get a list of all successfully built models from a.outputs
    ok_models = [x for x in a.outputs if x['failure'] is None]
    # Rank the models by DOPE score
    key = 'DOPE score'
    if sys.version_info[:2] == (2,3):
        # Python 2.3's sort doesn't have a 'key' argument
        ok_models.sort(lambda a,b: cmp(a[key], b[key]))
    else:
        ok_models.sort(key=lambda a: a[key])
    '''
    m = ok_models[0]

    print("Top model: %s (DOPE score %.3f)" % (m['name'], m[key]))

    # Get top 10 models
    for q in range(0,10):
    	m = ok_models[q]
    	print("Top model: %s (DOPE score %.3f)" % (m['name'], m[key]))
    '''

if __name__ == '__main__':
    main()
