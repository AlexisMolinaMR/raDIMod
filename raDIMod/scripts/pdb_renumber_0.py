def renumber_pdb(pdb_name, pdb_path, current_path):
    import os
    import sys
    import subprocess
    from SBI.structure.chain import Chain
    from SBI.sequence import Sequence
    from SBI.structure import PDB
    from Bio import SeqIO
    from Bio import ExPASy
    from Bio import AlignIO
    from Bio.Align import Applications

    clustal_exe = "/soft/EB_repo/bio/sequence/programs/goolf/1.7.20/ClustalW2/2.1/bin/clustalw2" #path where ClsutalW2 is stored

    name_pdb = ".".join(pdb_name.split('/')[-1].split('.')[:-1]) #retrieve of pdb code out of name

    new_pdb = PDB() #initialize a blank pdb file

    pdb_file = os.path.join(current_path,pdb_name) #full path where pdb file (template) is stored
    pdb = PDB(pdb_file) #initialize pdb routine from SBI library
    pdb.clean() # SBI's cleaning function

    for chain_id,chain_seq in sequences.iteritems():
        name_chain = current_path + "/" + name_pdb + "_" + chain_id
        name_seq   = current_path + "/" + chain_seq.get_identifier()
        pdb_chain  = current_path + "/" + pdb.get_chain_by_id(chain_id)
        new_chain  = current_path + "/" + Chain(name_pdb,chain_id)

        #define/create files
        infile = name_chain + "_" + name_seq + ".fa" #fasta file with the two sequneces for ClustalW2 to work with
        outfile = name_chain + "_" + name_seq + ".aln" #alignment file ot the two sequences
        dndfile = name_chain + "_" + name_seq + ".dnd" #tree outputted by ClustalW2 in Newick format

        #write the fasta infile
        fd = open(infile, "w") #open fasta file
        fd.write(">{0:s}\n{1:s}\n".format(name_chain,pdb_chain.protein_sequence)) #write both sequences in the input file
        fd.write(">{0:s}\n{1:s}\n".format(name_seq,chain_seq.get_sequence()))
        fd.close()

        try:
            #from now we use BioPython
            # run clustalw2
            msa_cline = Applications.ClustalwCommandline(clustal_exe, infile = infile, outfile = outfile)
            child = subprocess.Popen(str(msa_cline), stdout = subprocess.PIPE, stderr = subprocess.PIPE, shell = "/bin/bash")
            child.communicate()
            #store alignment in compare
            alignment = AlignIO.read(outfile,'clustal')
            structure = alignment[0].seq
            reference = alignment[1].seq
            try:
              len_3d = len(structure)
              len_ref = len(reference)
            except Exception as e:
              sys.stderr.write("ERROR: %s\n"%e)
              return e
        except Exception as e:
            sys.stderr.write("ERROR: %s\n"%e)
            return e
        #mapping of residues to the original sequence
        mapping = create_mapping(pdb_chain.protein_idx.split(";"),structure,reference)
        #fill the new chain with the correct numbering of residues
        for residue in pdb_chain.aminoacids:
           pair = (str(residue.number),residue.version)
           number,version = mapping.get(pair)
           residue.number = number
           residue.version = version
           new_chain.add_residue(residue)
        #fill the new pdb
        new_pdb.add_chain(new_chain)
    return new_pdb, alignment

print(renumber_pdb('2H5Y_B_139.pdb', '/home/amolina/Desktop/1ATG_A_build_trial', '/home/amolina/Desktop/1ATG_A_build_trial'))
