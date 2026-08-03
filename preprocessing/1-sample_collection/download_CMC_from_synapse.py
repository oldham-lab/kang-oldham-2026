import os
import synapseclient
import synapseutils
import argparse
import getpass

parser = argparse.ArgumentParser(description='Download RNAseq and genotypes of CMC.')
#parser.add_argument('RNAseq_directory', help='Directory to download RNAseq data to')
#parser.add_argument('Genotype_directory', help='Directory to download genotypes to')

args = parser.parse_args()


user = input("Synapse username:")
password = getpass.getpass('Synapse password:')

syn = synapseclient.Synapse()
syn.login(user,password)

print('Sync CMC')

# RNAseq release 3 fastq
files = synapseutils.syncFromSynapse(syn, 'syn18134196', path = 'RNAseq')
# Genotypes
#files = synapseutils.syncFromSynapse(syn, 'syn3275211', path = args.Genotype_directory)
# Metadata
files = synapseutils.syncFromSynapse(syn, 'syn3354385', path = 'metadata/') # Clinical metadata
files = synapseutils.syncFromSynapse(syn, 'syn3346807', path = 'metadata/') # Release 1 metadata
files = synapseutils.syncFromSynapse(syn, 'syn18358379', path = 'metadata/') # Release 3 metadata
files = synapseutils.syncFromSynapse(syn, 'syn24173489', path = 'metadata/') # Release 4 metadata
files = synapseutils.syncFromSynapse(syn, 'syn35641127', path = 'metadata/') # Release 6 metadata, ACC data only
