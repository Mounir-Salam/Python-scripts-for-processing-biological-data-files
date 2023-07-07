import os
import re
import shutil
import stat

# Set Paths
Directory = "C:\\Users\\user\\Documents\\Capstone_Project\\Conservation_Study"
Input = os.path.join(Directory, 'seqdump_updated.txt')
Output = os.path.join(Directory, "seqdump_duplicates")

# Create new output txt file
if os.path.exists(Output):
    os.chmod(Output, stat.S_IWUSR)
    shutil.rmtree(Output)
Output = open(os.path.join(Directory, "seqdump_duplicates.txt"), 'w+')

# Read the blast file to extract needed information
with open(Input, 'r+') as file:
    line = file.readline()
    UniprotID = []
    info = []
    dictionary = {}
    while(line != ''):
        if line.startswith('>'):
            replace = '>'
            
            # Extract protein name and organism name
            names = re.search('(\{.*?\})$', line).group(1)

            #extract Uiprot ID
            code = re.search('>(.*?)\{', line).group(1)

            if names in info:
                if names in dictionary.keys():
                    dictionary[names].append(code)
                else:
                    dictionary[names] = [UniprotID[info.index(names)], code]
            else:
                info.append(names)
                UniprotID.append(code)

        line = file.readline()

    for key in dictionary.keys():
        Output.write(key + ": " + str(dictionary[key]) + "\n")

Output.close()