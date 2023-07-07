import os
import re
import shutil
import stat

# set paths
PDB_Directory = "C:/Users/user/Documents/Capstone_Project/MAPK14_PDB_IDs/Output2"
Index_Directory = "C:/Users/user/Documents/Capstone_Project/PDBbind_v2020_plain_text_index.tar/index"

Index_Input = os.path.join(Index_Directory, "INDEX_general_PL.2020")
Output = os.path.join(PDB_Directory, "Affinity_Data")

# Create new output csv file
if os.path.exists(Output):
    os.chmod(Output, stat.S_IWUSR)
    shutil.rmtree(Output)
Output = open(os.path.join(PDB_Directory, 'Affinity_Data.csv'), 'w+')
Output.write('PDB ID,Resolution,Ligand,Affinity\n')

# Extract PDB IDs from the directory of your concern
PDB_list = []
for pdb in os.listdir(PDB_Directory):
    if (pdb.endswith('.ent')):
        PDB_list.append(re.search('pdb(.+?)\.ent', pdb).group(1))

# Read Index File and extract needed affinity data
with open(Index_Input, 'r+') as file:
    line = file.readline()
    while(line != ''):
        if line.startswith('#'):
            line = file.readline()
            continue
        else:
            line = line.strip().split(maxsplit=6)
            if line[0] in PDB_list:
                Output.write(line[0] + ',' + line[1] + ',' + re.search('\((.+?)\)', line[6]).group(1) + ',' + line[3] + '\n')
            line = file.readline()

Output.close()
            