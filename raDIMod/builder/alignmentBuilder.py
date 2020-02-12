import argparse as ap
import os
import sys

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
    file = open(path).readlines()
    print("Opening " + path + '\n')
    for lines in file:
        k = lines.split()
        if "#" not in k[0][0] and k[0] not in codes_arc:
            print("Retrieving code: " + k[0])
            codes_arc.append(k[0])

    print("{} codes retrieved!\n".format(len(codes_arc)))
    return codes_arc

def retrieve_sequence(path,codes_arc):
    file = open(path).readlines()
    print("Retrieving sequence...")
    for lines in file:
        k = lines.split()
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

def retrieve_template_codes(path):
    file = open(path).readlines()
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

    return temp_codes


def aliBuild(codes_arc, seq, temp_codes):

    print("Building alignment...")

    seq += "*"
    f1 =  open("ali.pir", "w")

    f1.write(">P1;{}\n".format(codes_arc[0]))
    f1.write("SEQUENCE:{}: : : : : : : SEQUENCE :(SEQU) .\n".format(codes_arc[0]))
    f1.write(seq+'\n')
    for i in temp_codes:
        f1.write(">P1;{}\n".format(i[0]))
        f1.write("structureX:{}:{}:{}:{}:::::\n".format(i[0],i[2][0],i[2][2],i[2][1]))
        f1.write(i[1]+'\n')

    print("ali.pir built successfully!\n")

def main():

    path = parseArg()

    retrieve_codes(path=path)
    retrieve_sequence(path=path,codes_arc=codes_arc)
    retrieve_template_codes(path=path)
    retrieve_relative_postions(temp_codes=temp_codes)
    fill_temp_codes_tuple(temp_codes=temp_codes, rel_pos=rel_pos)
    aliBuild(codes_arc=codes_arc,seq=retrieve_sequence(path=path,codes_arc=codes_arc),temp_codes=temp_codes)

    print("Ahead with your models!")


if __name__ == "__main__":
    main()
