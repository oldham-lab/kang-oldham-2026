import os
import synapseclient
import synapseutils
import argparse
import getpass 

#parser = argparse.ArgumentParser(description='Download RNAseq and genotypes of CMC.')
#parser.add_argument('RNAseq_directory', help='Directory to download RNAseq data to')

#args = parser.parse_args()

user = input("Synapse username:")
password = getpass.getpass('Synapse password:')

syn = synapseclient.Synapse()
syn.login(user,password)

print('Sync Brainseq')
# RNAseq
files = synapseutils.syncFromSynapse(syn, 'syn15672826', path = 'RNAseq/', ifcollision = 'keep.local')

# Phenotype file
files = synapseutils.syncFromSynapse(syn, 'syn4566330', path='phenotype_data/')
#files = synapseutils.syncFromSynapse(syn, 'syn7203084', path='phenotype_data/') #LIBD_WGBS whole genome bisulfite
#files = synapseutils.syncFromSynapse(syn, 'syn7203089', path='phenotype_data/') # snpArray_LIBD-WGBS_metadata.csv
#files = synapseutils.syncFromSynapse(syn, 'syn8017780', path='phenotype_data/') # bisulfiteSeq_LIBD-WGBS_metadata.csv
