import pandas as pd
import numpy as np
import os
import shutil
import stat
import xml.etree.cElementTree as et

# Specify directories containing needed data
Directory_XML = 'C:/Users/user/Documents/Capstone_Project/Protein-Ligand_Data_PLIP_XML'
Directory_Ligands = 'C:/Users/user/Documents/Capstone_Project/MAPK14_PDB_IDs/Output'

# Extracting needed ligands
Ligands = pd.read_csv(os.path.join(Directory_Ligands, 'Interaction_info.csv'))
Ligands = pd.array(Ligands.get('Ligands'))

#####################################################################
# Create new output csv files

## Output files to be created
Output_Hydrophobic = os.path.join(Directory_XML, 'Interaction_info_Hydrophobic_Interactions')
Output_HBond = os.path.join(Directory_XML, 'Interaction_info__Hydrogen_Bonds')
Output_piStack = os.path.join(Directory_XML, 'Interaction_info_pi_Stacking')

## For Hydrophobic Interactions
if os.path.exists(Output_Hydrophobic):
    os.chmod(Output_Hydrophobic, stat.S_IWUSR)
    shutil.rmtree(Output_Hydrophobic)
# Output_Hydrophobic = open(os.path.join(Directory_XML, 'Interaction_info_Hydrophobic_Interactions.csv'), 'w+')
# Output_Hydrophobic.write('PDB ID,Ligand,Residue,Distance,Ligand Atom,Protein Atom\n')
Hydrophobic_df = pd.DataFrame(columns=['PDB ID','Ligand','Residue','Distance','Ligand Atom','Protein Atom'])

## For Hydrogen Bonds
if os.path.exists(Output_HBond):
    os.chmod(Output_HBond, stat.S_IWUSR)
    shutil.rmtree(Output_HBond)
# Output_HBond = open(os.path.join(Directory_XML, 'Interaction_info__Hydrogen_Bonds.csv'), 'w+')
# Output_HBond.write('PDB ID,Ligand,Residue,Distance (D-A),Protein Donor?,Side Chain?,Donor Atom,Acceptor Atom\n')
HBond_df = pd.DataFrame(columns=['PDB ID','Ligand','Residue','Distance (D-A)','Protein Donor?','Side Chain?','Donor Atom','Acceptor Atom'])

## For pi-Stacking
if os.path.exists(Output_piStack):
    os.chmod(Output_piStack, stat.S_IWUSR)
    shutil.rmtree(Output_piStack)
# Output_piStack = open(os.path.join(Directory_XML, 'Interaction_info_pi_Stacking.csv'), 'w+')
# Output_piStack.write('PDB ID,Ligand,Residue,Distance,Stacking Type,Ligand Atoms\n')
piStack_df = pd.DataFrame(columns=['PDB ID','Ligand','Residue','Distance','Stacking Type','Protein Atoms','Ligand Atoms'])
#####################################################################

# Read xml files
for xml in os.listdir(Directory_XML):
    if (xml.endswith('.xml')):
        Ligands_tree = et.parse(os.path.join(Directory_XML, xml))
        root = Ligands_tree.getroot()
        PDB = root.find('pdbid').text
        print(PDB)
        
        for binding_Site in root.iter('bindingsite'):
            Ligand = binding_Site.find('identifiers/longname').text
            if Ligand not in Ligands:
                continue
            else:                
                # Hydrophobic Interactions
                for interaction in binding_Site.iter('hydrophobic_interaction'):
                    # PDB ID = <pdbid>
                    # Ligand = <longname>
                    # Residue = <restype> + <resnr> + ' ' + <reschain>
                    # Distance = <dist>
                    # Ligand Atom = <ligcarbonidx>
                    # Protein Atom = <protcarbonidx>
                    
                    dfFind = interaction.find

                    info = [PDB, Ligand,
                            dfFind('restype').text + dfFind('resnr').text + ' ' + dfFind('reschain').text,
                            dfFind('dist').text, dfFind('ligcarbonidx').text,
                            dfFind('protcarbonidx').text]
                    Hydrophobic_df.loc[len(Hydrophobic_df)] = info
                
                # Hydrogen Bonds
                for interaction in binding_Site.iter('hydrogen_bond'):
                    # PDB ID = <pdbid>
                    # Ligand = <longname>
                    # Residue = <restype> + <resnr> + ' ' + <reschain>
                    # Distance (D-A) = <dist_d-a>
                    # Protein Donor? = <protisdon>
                    # Side Chain? = <sidechain>
                    # Donor Atom = <donortype>
                    # Acceptor Atom = <acceptortype>

                    dfFind = interaction.find

                    info = [PDB, Ligand,
                            dfFind('restype').text + dfFind('resnr').text + ' ' + dfFind('reschain').text,
                            dfFind('dist_d-a').text, dfFind('protisdon').text,
                            dfFind('sidechain').text, dfFind('donortype').text,
                            dfFind('acceptortype').text]
                    HBond_df.loc[len(HBond_df)] = info

                # pi-Stackings
                for interaction in binding_Site.iter('pi_stack'):
                    # PDB ID = <pdbid>
                    # Ligand = <longname>
                    # Residue = <restype> + <resnr> + ' ' + <reschain>
                    # Distance = <centdist>
                    # Stacking Type = <type>
                    # Protein Atoms = <idx> + ',' + <idx> + ... (in <prot_idx_list>)
                    # Ligand Atoms = <idx> + ',' + <idx> + ... (in <lig_idx_list>)

                    dfFind = interaction.find
                    Protein_atoms = []
                    Ligand_atoms = []

                    # Extract Protein Atoms
                    for atom in dfFind('prot_idx_list').iter('idx'):
                        Protein_atoms.append(atom.text)
                    # Extract Ligand Atoms
                    for atom in dfFind('lig_idx_list').iter('idx'):
                        Ligand_atoms.append(atom.text)

                    info = [PDB, Ligand,
                            dfFind('restype').text + dfFind('resnr').text + ' ' + dfFind('reschain').text,
                            dfFind('centdist').text, dfFind('type').text,
                            ', '.join(Protein_atoms), ', '.join(Ligand_atoms)]
                    piStack_df.loc[len(piStack_df)] = info

# Create Output .csv files
Hydrophobic_df.to_csv(Output_Hydrophobic + '.csv')
HBond_df.to_csv(Output_HBond + '.csv')
piStack_df.to_csv(Output_piStack + '.csv')
