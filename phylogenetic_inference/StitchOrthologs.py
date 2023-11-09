#!/bin/env python3
import sys
import re
import pandas as pd
import numpy as np
from os.path import exists
from io import StringIO
from Bio import AlignIO
from Bio import SeqIO
from Bio import Phylo

MAF=sys.argv[1]
OUTPUT_FASTA=sys.argv[2]
OUTPUT_CSV=sys.argv[3]

if exists(OUTPUT_FASTA): #Write MAF depth stats to a file if that file does not already exist (this should be changed with a command line argument named flag option)
	print("\n\t"+OUTPUT_FASTA+" exists. Please rename and try again.")
	print("\n\tUSAGE: StitchOrthologs.py [aln.maf] [out.fasta] [out.csv]\n")
	exit()

if exists(OUTPUT_CSV): #Write MAF depth stats to a file if that file does not already exist (this should be changed with a command line argument named flag option)
	print("\n\t"+OUTPUT_CSV+" exists. Please rename and try again.")
	print("\n\tUSAGE: StitchOrthologs.py [aln.maf] [out.fasta] [out.csv]\n")
	exit()

#########################################################################################
# FUNCTION: IS BLOCK A 1:1 ORTHOLOG?
#########################################################################################

def check_block_orthology(MAF_BLOCK,TAXA_LIST): #[ARG1=AlignIO MAF Block], [ARG2=List of tip taxa in MAF tree topology]
	SPECIES_IN_BLOCK=[]
	if len(MAF_BLOCK) != len(TAXA_LIST):	
		return False
	else:
		for SEQ in MAF_BLOCK: 
			if len(SEQ.id.split(".")) == 2: #For [species].[seq] format
				if SEQ.id.split(".",2)[0] not in SPECIES_IN_BLOCK: #If seqrecord is the first occurance of the species in this block, write their name to the SPECIES_IN_BLOCK list
					SPECIES_IN_BLOCK.append(SEQ.id.split(".",2)[0])
			else:
				if ".".join(SEQ.id.split(".",2)[0:2]) not in SPECIES_IN_BLOCK: #If seqrecord is the first occurance of the species in this block, write their name to the SPECIES_IN_BLOCK list
					SPECIES_IN_BLOCK.append(".".join(SEQ.id.split(".",2)[0:2]))
		
		# print("Species in block")
		# print(sorted(SPECIES_IN_BLOCK))
		# print("Taxa list")
		# print(sorted(TAXA_LIST))
		
		if sorted(SPECIES_IN_BLOCK) == sorted(TAXA_LIST):
			return True
		else:
			return False

#########################################################################################
# READ TREE TOPOLOGY AND LOAD TIPS INTO A LIST
#########################################################################################

with open(MAF) as file: #Tree topology is stored as second line in MAF
	next(file)
	TREE=next(file)

TREE=StringIO(re.sub(r"^.* \(","(",TREE).strip()) #Remove whitespace and parse string with StringIO so it can be passed to Phylo.read()
TREE=Phylo.read(TREE,"newick")
LEAVES=sorted([TIP.name for TIP in TREE.get_terminals()]) #Load tree tips into a list for generating aln metrics

#########################################################################################
# SETUP DATAFRAME COLUMNS
#########################################################################################

DF_COLS=["BlockID","AlnLength","AlnStart","AlnStop"]
for TAXON in LEAVES:
	DF_COLS.append(str(TAXON+"Start"))
	DF_COLS.append(str(TAXON+"Stop"))
	DF_COLS.append(str(TAXON+"Polarity"))
INDICES_LIST=[]

#########################################################################################
# GENERATE 1:1 ORTHOLOGY INDEX
#########################################################################################

BLOCK_NUM=0
ALN_START=0
for BLOCK in AlignIO.parse(MAF,"maf"):
	if BLOCK_NUM == 0:
		REF_SEQ=BLOCK[0].id
		REF_GENOME=BLOCK[0].id.split(".",2)[0] #Pull reference genome from first block (first SeqRecord in first block)
	if check_block_orthology(BLOCK,LEAVES) == True: #Skips any alignment blocks evaluated to be anything other than a 1:1 ortholog block
		#CREATES INDEX LINE
		ROW_DICT={} #Dictionary to be populated to transfer to ALN_INDEX DataFrame
		ROW_DICT['BlockID']=BLOCK_NUM #Populate dictionary with BlockID
		ROW_DICT['AlnLength']=len(BLOCK[0].seq) #Populate dictionary with AlnLength
		ROW_DICT['AlnStart']=ALN_START #Populate dictionary with AlnStart
		ROW_DICT['AlnStop']=ALN_START+len(BLOCK[0].seq)-1
		
		for SEQREC in BLOCK: #Loop through sequence records in BLOCK
			if len(SEQREC.id.split(".")) == 2: #For [species].[seq] format
				TAXON=SEQREC.id.split(".",2)[0]
			else:
				TAXON=".".join(SEQREC.id.split(".",2)[0:2])
			START=SEQREC.annotations['start']
			POLARITY=SEQREC.annotations['strand']
			if POLARITY == -1: #If sequence polarity is -1, STOP=START-Length
				STOP=START-SEQREC.annotations['size']+1
			else: #If sequence polarity is 1, STOP=START+Length
				STOP=START+SEQREC.annotations['size']-1
			ROW_DICT[str(TAXON+"Start")]=START
			ROW_DICT[str(TAXON+"Stop")]=STOP
			ROW_DICT[str(TAXON+"Polarity")]=POLARITY

		ROW_LIST=[[ROW_DICT[KEY] for KEY in DF_COLS]]		
		INDICES_LIST+=ROW_LIST

		ALN_START+=(len(BLOCK[0].seq)) #Modify ALN_START for next block

	BLOCK_NUM+=1

ALN_INDEX=pd.DataFrame(INDICES_LIST,columns=DF_COLS) #Create pandas dataframe from 1:1 ortholog block data
ALN_INDEX.to_csv(OUTPUT_CSV,index=False)
print("Alignment index CSV written to: "+OUTPUT_CSV)

#########################################################################################
# GENERATE SPLICED ALIGNMENT
#########################################################################################

MAF_INDEX=AlignIO.MafIO.MafIndex(sqlite_file=(MAF+".sqlite3.idx"),maf_file=MAF,target_seqname=REF_SEQ) #Create index for MAF
STARTS=ALN_INDEX[REF_GENOME+"Start"].tolist()
STOPS=[x+1 for x in ALN_INDEX[REF_GENOME+"Stop"].tolist()] #get_spliced function uses an exclusive stop index, whereas the list generated above used the standard python list syntax of inclusive end coordinates

SPLICED_ALN=MAF_INDEX.get_spliced(STARTS,STOPS) #Generate spliced alignment; rate-limiting step
SeqIO.write(SPLICED_ALN,OUTPUT_FASTA,"fasta")

#########################################################################################
# STD_OUT PRINTING
#########################################################################################

print("Alignment matrix FASTA written to: "+OUTPUT_FASTA)
print("\tMatrix dimensions: %i x %i" % (SPLICED_ALN.get_alignment_length(),len(SPLICED_ALN),))