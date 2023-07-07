import os
import re
import shutil
import stat

# Set Paths
Directory = "C:\\Users\\user\\Documents\\Capstone_Project\\Conservation_Study"
Input = os.path.join(Directory, 'seqdump.txt')
Output = os.path.join(Directory, "seqdump_updated")

# Create new output txt file
if os.path.exists(Output):
    os.chmod(Output, stat.S_IWUSR)
    shutil.rmtree(Output)
Output = open(os.path.join(Directory, "seqdump_updated.txt"), 'w+')

# Read the blast file to extract needed information
with open(Input, 'r+') as file:
    line = file.readline()
    while(line != ''):
        if line.startswith('>'):
            replace = '>'
            
            #extract Uiprot ID
            search = re.search('>(.*?)\.', line)
            replace = replace + search.group(1)

            # Extract short name of the protein
            # Extract full name if short name unavailable
            search = re.search('^.*?\\bShort=\\b(.*?)[;\[]', line)
            if search == None:
                search = re.search('^.*?\\bFull=\\b(.*?)[;\[]', line)
            
            # For better format, we avoid spaces since some webtools takes spaces as a stopping criteria
            # Replace spaces with underscores   
            search = re.sub(' ', '_', search.group(1).strip())
            replace = (replace + '{'+ search + '}')

            # Extract organism name (Taxonomy)
            search = re.search('\[(\w)\w+ (\w+).*\]', line)
            search = search.group(1)+"."+search.group(2)
            replace = replace+'{'+search+'}\n'

            # The components of 'replace' entity is appended to the new output file instead of the old line
            Output.write(replace)
        
        # Lines that we dont need to extract any information from
        # The 'line' entity would in this case contain the protein sequence
        # Will be appended directly without any alterations
        else:
            Output.write(line)

        line = file.readline()

Output.close()

