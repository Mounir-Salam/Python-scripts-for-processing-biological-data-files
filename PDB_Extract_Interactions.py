import os
import shutil
import stat
import re

# Set paths
Directory = "C:/Users/user/Documents/Capstone_Project/MAPK14_PDB_IDs/Output2"
Output = os.path.join(Directory, "Interaction_info")

# Create new output csv file
if os.path.exists(Output):
    os.chmod(Output, stat.S_IWUSR)
    shutil.rmtree(Output)
Output = open(os.path.join(Directory, "Interaction_info.csv"), 'w+')
Output.write("PDB ID,Ligands,Interaction Sites\n")

# Ligands we can ignore
ignore = ['SO4','GOL','MG','ACT','CL','ZN','NA','CSO','CSS','CME','TPO','PTR','DMS','BOG','PO2','EDO','BME','PEG']
pocket = ['VAL30 A','GLY31 A','SER32 A','TYR35 A','VAL38 A','ALA51 A','VAL52 A','LYS53 A','ARG67 A','ARG70 A','GLU71 A','LEU74 A','LEU75 A','VAL83 A','ILE84 A','LEU86 A','LEU104 A','VAL105 A','THR106 A','HIS107 A','LEU108 A','MET109 A','GLY110 A','ALA111 A','ASP112 A','ILE141 A','ILE146 A','HIS148 A','SER154 A','ALA157 A','LEU167 A','ASP168 A','PHE169 A','LEU171 A']

for pdb in os.listdir(Directory):
    if (pdb.endswith(".ent")):
        pdbPath = os.path.join(Directory,pdb)
        #pdbcsv = re.search('^pdb(.*?)\.', pdb)
        with open(pdbPath, 'r+') as filePDB:

            line = filePDB.readline()
            ligands = {}
            interactions = {}

            while (line != ""):
                if line.startswith('REMARK 800'):
                    #print(pdb)
                    
                    # Correlate active site code with its respective ligands
                    if 'SITE_IDENTIFIER' in line:

                        # extract active site code (AC1, AC2, ...)
                        line = line.split()
                        identifier = line[-1]

                        # navigate to 'SITE_DESCRIPTION' 2 lines after 'SITE_IDENTIFIER'
                        # extract name of ligand
                        filePDB.readline()
                        line = filePDB.readline()
                        mol = re.search('RESIDUE (.+?) ', line, re.IGNORECASE).group(1)

                        # Store in ligands dictionary format {identifier:mol}
                        if mol not in ignore:
                            ligands[identifier] = mol
                        
                        line = filePDB.readline()

                if line.startswith('SITE'):
                    line = line.split()
                    if line[2] in ligands.keys():
                        # Interactions begin at index 4
                        i=4
                        amino_acid = ''
                        chain = ''
                        position = ''

                        # Create new entry in dictionary for each valid ligand
                        if line[2] not in interactions.keys():
                            interactions[line[2]] = []

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
                            interactions[line[2]].append(amino_acid + str(position) + " " +chain)

                line = filePDB.readline()

            # Check if the ligand is in the ATP-binding pocket
            # And append to ouptut 
            for key in ligands.keys():
                for residue in interactions[key]:
                    if residue in pocket:
                        #Output.write(pdbcsv.group(1) + ',' + ligands[key] + ",\"" + ", ".join(interactions[key]) + "\"\n")
                        Output.write(pdb + ',' + ligands[key] + ",\"" + ", ".join(interactions[key]) + "\"\n")
                        break
        
Output.close()