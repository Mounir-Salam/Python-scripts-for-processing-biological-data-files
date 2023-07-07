import os
import shutil
import stat
import re

Directory = "C:/Users/user/Documents/Capstone_Project/MAPK14_PDB_IDs"

def filterPDB(Directory):

    def normalize(Length, Position):
        return (0.03-0.02)*(abs(Position-Length/2)/(Length-length/2)) + 0.02

    count = 0
    resolution = 0

    # Create Directory 'Output' to Store Valid PDB Structures
    Output = os.path.join(Directory, 'Output')
    if os.path.exists(Output):
        os.chmod(Output, stat.S_IWUSR)
        shutil.rmtree(Output)
    os.mkdir(Output)
    os.chmod(Output, stat.S_IWUSR)
    info = open(os.path.join(Output, "info.csv"), 'w+')
    info.write('PDB File,Resolution,Mutation\n')

    for pdb in os.listdir(Directory):
        if (pdb.endswith(".ent")):
            pdbPath = os.path.join(Directory,pdb)
            mutant = 'NO'
            unfit = False
            with open(os.path.join(Directory, pdb), 'r+') as filePDB:

                line = filePDB.readline()
                while (not unfit) and (line != ""):
                    
                    # Check if Structure is a mutant on the original
                    if line.startswith('COMPND') and 'MUTATION:' in line:
                        mutant = re.sub('[^a-zA-Z]','',line.split()[3])

                    # Check if Experimental Method is X-Ray 
                    if line.startswith('EXPDTA'):
                        if "X-RAY DIFFRACTION" not in line:
                            unfit = True
                            break

                    # Check if Resolution (i.e REMARK 2) <= 2.5 Angstroms
                    if line.startswith('REMARK   2'):
                        line = (filePDB.readline()).split()
                        if (eval(line[3]) <= 2.5):
                            resolution = line[3]
                            line = filePDB.readline()
                        else: 
                            unfit = True
                            break

                    # Collect Missing Residues data (i.e REMARK 465)
                    if line.startswith('REMARK 465'):
                        value = 465
                        missingResidues = {}

                        for i in range(6):
                            filePDB.readline()

                        line = (filePDB.readline()).split()
                        value = eval(line[1])

                        while (value == 465):
                            chain = line[3]

                            if chain not in missingResidues:
                                missingResidues[chain] = []

                            if eval(line[4]) > 0:
                                missingResidues[chain].append(eval(line[4]))

                            line = (filePDB.readline()).split()
                            value = eval(line[1])

                        line = filePDB.readline()

                    # Find length of polypeptide chain
                    # Filter out structures with alot of missing residues
                    if line.startswith('DBREF') and line.split()[2] in missingResidues:
                        print(pdb)
                        line = line.split()
                        length = eval(line[4])
                        chain = line[2]

                        #print(pdb, chain, missingResidues)
                        if len(missingResidues[chain]) > length*0.1:
                            #print('unfit1')
                            unfit = True
                            break

                        if len(missingResidues[chain]) > 0:
                            #Position = 0
                            gap = 0
                            start = 0
                            end = 0

                            for residue in missingResidues[chain]:
                                # use normalization formula to turn residue distance from the center into a percentage
                                # of allowed missing residues between 0.02 and 0.03 such that residues farthest from the
                                # center are allowed to have 3% missing residues from the entire length of the polypeptide
                                # and residues at the center are allowed 2% becuase missing residues in the center are more
                                # probable to have functional importance when compared to residues at the C and N terminus
                                if start == 0:
                                    start = end = residue
                                
                                # if the gap is done or if we reach the end of the sequence
                                # check length of gap
                                # check if missing residues are only at the beginning and end of sequence
                                if residue != end+1 or residue == missingResidues[chain][len(missingResidues[chain])-1]:
                                    #Position = (start+end)/2
                                    gap = end-start
                                    # print(Position, round(length*normalize(length, Position)), start, end, gap)

                                    # Normalization Formula
                                    #if gap > round(length*normalize(length, Position)):
                                        #unfit = True
                                        #break
                                    if (start == 1) or (residue == length):
                                        if gap > length*0.05:
                                            #print('unfit2')
                                            unfit = True
                                            break
                                        
                                        else:
                                            start = end = residue
                                    
                                    else:
                                        #print('unfit3', start, end, gap)
                                        unfit = True
                                        break
                                            
                                else: end = residue

                            if unfit: break

                    line = filePDB.readline()

                # Copy PDB Structure to Output Directory
                if not unfit:
                    count += 1
                    info.write(pdb + "," + resolution + "," + mutant + '\n')
                    shutil.copy(pdbPath, Output)
                
    print(count)
    info.close()
 
filterPDB(Directory)