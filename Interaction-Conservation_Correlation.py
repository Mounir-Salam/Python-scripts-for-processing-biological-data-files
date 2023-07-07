import os
import re
import shutil
import stat

Directory1 = 'C:/Users/user/Documents/Capstone_Project/MAPK14_PDB_IDs/Output2'
Directory2 = 'C:/Users/user/Documents/Capstone_Project/Conservation_Study/Consurf_Results/Consurf_Outputs_5XYY_Final'

Input1 = os.path.join(Directory1, 'Residues_Stat_no_extra.csv')
Input2 = os.path.join(Directory2, 'msa_aa_variety_percentage.csv')
Output = os.path.join(Directory1, 'Correlation')

# Create new output csv file
if os.path.exists(Output):
    os.chmod(Output, stat.S_IWUSR)
    shutil.rmtree(Output)
Output = open(os.path.join(Directory1, 'Correlation.csv'), 'w+')
Output.write('Residue,Number of Occurrences,MAX AA,ConSurf Grade\n')

with open(Input1, 'r+') as stats, open(Input2, 'r+') as conservation:
    for i in range(4):
        conservation.readline()

    # Extract column names    
    scolumns = stats.readline().strip().split(',')
    ccolumns = conservation.readline().strip().split(',')
    
    # Extract indecies of needed columns
    residue_index = scolumns.index('Residue')
    pos_index = ccolumns.index('pos')
    aa_index = ccolumns.index('MAX AA')
    score_index = ccolumns.index('ConSurf Grade')

    sline = stats.readline().strip()
    cline = conservation.readline().strip()
    cline = cline.split(',')

    while sline != '':
        sline = sline.split(',')
        # Extract residue position from input1
        residue = re.search('\d+', sline[residue_index]).group(0)
        while cline != '':
            # Extract residue position from input2
            pos = cline[pos_index]

            # In the end the Output file will only contain the residues in the position that contributed in the interactions
            if residue == pos:
                Output.write(','.join(sline) + ',' + cline[aa_index] + ',' + cline[score_index] + '\n')
                #cline = conservation.readline().strip()
                break
            cline = conservation.readline().strip()
            cline = cline.split(',')
        sline = stats.readline().strip()
        
        # Ignore water interaction from input1
        if 'HOH' in sline:
            sline = stats.readline().strip()

Output.close()