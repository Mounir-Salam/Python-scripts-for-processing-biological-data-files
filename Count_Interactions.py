import os
import re
import shutil
import stat

# Set Paths
Directory = "C:/Users/user/Documents/Capstone_Project/MAPK14_PDB_IDs/Output2"
Input = os.path.join(Directory, 'Interaction_info.csv')
Output = os.path.join(Directory, 'Residues_Stat')

# Create new output csv file
if os.path.exists(Output):
    os.chmod(Output, stat.S_IWUSR)
    shutil.rmtree(Output)
Output = open(os.path.join(Directory, 'Residues_Stat.csv'), 'w+')
Output.write('Residue,Number of Occurrences\n')

# Read the blast file to extract needed information
with open(Input, 'r+') as file:
    file.readline()
    line = file.readline()
    counter = {}
    while (line != ''):
        line = re.sub('\"|(?<=,) |\n', '', line)
        line = line.split(',', 2)[2].split(',')
        for residue in line:

            if 'HOH' in residue: 
                residue = 'HOH'
            
            if residue in counter.keys():
                counter[residue] += 1
            else:
                counter[residue] = 1

        line = file.readline()

    # Sort list based on residue position
    keys = list(counter.keys())
    if 'HOH' in keys:
        keys.remove('HOH')
        keys.sort(key= lambda x: int(re.search('\d+', x).group(0)))
        keys.append('HOH')
    else:
        keys = keys.sort(key= lambda x: int(re.search('\d+', x).group(0)))
    
    for key in keys:
        Output.write(key + ',' + str(counter[key]) + '\n')

Output.close()