import os
import shutil
import stat
import re

# set paths
Directory = "C:/Users/user/Documents/Capstone_Project/MAPK14_PDB_IDs"
Output = os.path.join(Directory, 'Output2')

# Create Directory 'Output' to Store Valid PDB Structures
if os.path.exists(Output):
    os.chmod(Output, stat.S_IWUSR)
    shutil.rmtree(Output)
os.mkdir(Output)
os.chmod(Output, stat.S_IWUSR)
info = open(os.path.join(Output, "info.csv"), 'w+')
info.write('PDB File,Resolution,Mutation\n')

count = 0
resolution = 0
# Ligands we can ignore
ignore = ['SO4','GOL','MG','ACT','CL','ZN','NA','CSO','CSS','CME','TPO','PTR','DMS','BOG','PO2','EDO','BME','PEG']
# Residues of the pocket we want
pocket = ['VAL30 A','GLY31 A','SER32 A','TYR35 A','VAL38 A','ALA51 A','VAL52 A','LYS53 A','ARG67 A','ARG70 A','GLU71 A','LEU74 A','LEU75 A','VAL83 A','ILE84 A','LEU86 A','LEU104 A','VAL105 A','THR106 A','HIS107 A','LEU108 A','MET109 A','GLY110 A','ALA111 A','ASP112 A','ILE141 A','ILE146 A','HIS148 A','SER154 A','ALA157 A','LEU167 A','ASP168 A','PHE169 A','LEU171 A']

for pdb in os.listdir(Directory):
    if (pdb.endswith(".ent")):
        pdbPath = os.path.join(Directory,pdb)
        unfit = False
        ATP_Bind = False
        mutant = 'NO'
        with open(os.path.join(Directory, pdb), 'r+') as filePDB:

            line = filePDB.readline()
            ligands = {}

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

                # Collect Ligand Data
                if line.startswith('REMARK 800'):
                    while line.startswith('REMARK 800'):
                        # Correlate active site code with its respective ligands
                        if 'SITE_IDENTIFIER' in line:

                            # extract active site code (AC1, AC2, ...)
                            line = line.split()
                            identifier = line[-1]

                            # navigate to 'SITE_DESCRIPTION' 2 lines after 'SITE_IDENTIFIER'
                            # extract name of ligand
                            filePDB.readline()
                            line = filePDB.readline().split()
                            mol = line[-3]

                            # Store in ligands dictionary format {identifier:mol}
                            if mol not in ignore:
                                ligands[identifier] = mol
                            
                        line = filePDB.readline()
                
                # Find length of polypeptide chain
                # Filter out structures with alot of missing residues
                if line.startswith('DBREF') and line.split()[2] in missingResidues:

                    # DBREF is written somewhere after REMARK 800; check if ligands exist, if not -> unfit
                    if not ligands:
                        unfit = True
                        break

                    #print(pdb)
                    line = line.split()
                    length = eval(line[4])
                    chain = line[2]

                    #print(pdb, chain, missingResidues)
                    if len(missingResidues[chain]) > length*0.1:
                        #print('unfit1')
                        unfit = True
                        break
                    line = filePDB.readline()
                
                if line.startswith('SITE'):
                    line = line.split()
                    if line[2] in ligands.keys():
                        # Interactions begin at index 4
                        i=4
                        amino_acid = ''
                        chain = ''
                        position = ''

                        while i < len(line):
                            
                            # Extract amino acid name
                            if line[i].isalpha() and len(line[i])>1:
                                amino_acid = line[i]
                                i += 1
                            
                            # Extract chain written between amino acid and position
                            if line[i].isalpha() and len(line[i])==1:
                                chain = line[i]
                                i += 1

                            # Extract Position
                            if line[i].isnumeric():
                                position = line[i]
                                i += 1
                            
                            # Special case if the position is 4 digits it won't be separated from the chain by a space
                            else:
                                chain = line[i][0]
                                position = line[i][1:]
                                i += 1

                            residue = amino_acid + str(position) + " " +chain
                            if residue in pocket:
                                ATP_Bind = True
                                break


                line = filePDB.readline()

            # Copy PDB Structure to Output Directory
            if not unfit and ATP_Bind:
                count += 1
                info.write(pdb + "," + resolution + "," + mutant + '\n')
                shutil.copy(pdbPath, Output)
            
print(count)
info.close()