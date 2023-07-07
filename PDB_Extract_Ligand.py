import os
import shutil
import stat
import re

ligands = []

# Set paths
Directory = "C:/Users/user/Documents/Capstone_Project/MAPK14_PDB_IDs/Output2"
Output = os.path.join(Directory, "Ligand_info")

# Create new output csv file
if os.path.exists(Output):
    os.chmod(Output, stat.S_IWUSR)
    shutil.rmtree(Output)
Output = open(os.path.join(Directory, "Ligand_info.csv"), 'w+')
Output.write("PDB ID,Ligands\n")

for pdb in os.listdir(Directory):
    if (pdb.endswith(".ent")):
        pdbPath = os.path.join(Directory,pdb)
        pdbcsv = re.search('^pdb(.*?)\.', pdb)
        with open(os.path.join(Directory, pdb), 'r+') as filePDB:

            line = filePDB.readline()
            ligands = []
            while (line != ""):
                
                if line.startswith('HET  '):
                    #print(pdb)
                    line = line.split()
                    ligands.append(line[1])

                line = filePDB.readline()

            # Append to ouptut file
            #Output.write(pdbcsv.group(1) + ",\"" + ", ".join(ligands) + "\"\n")
            Output.write(pdb + ",\"" + ", ".join(ligands) + "\"\n")
        
Output.close()