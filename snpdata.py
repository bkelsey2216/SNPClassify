"""
	Jessica Jowdy, Brooke Kelsey
	CS68 Final Project, Fall 2014

	Program to run the first experimental method.
"""

from SNPClass import *
from DTree import *
import os

# Global Variables
# Dictionary storing all SNP objects in the chromosome range. Stores the population data.
SNPDict = {}
# List of all Patients patient objects. Stores each patient's SNP data within the chromosome range.
Patients = []
# List of all Node/Leaf objects added to the tree during the ID3 algorithm.
nodes = []


"""
	readPop() opens the population VCF file and reads the necessary information to create the individual SNP objects.
	The population VCF file is obtained by selecting "Download Genotype Data" -> "Data: Aggregates by Population" -> "Filter by: None or Populations with Checked Samples"
	It then calls makeSNPClass() to populate the SNPDict. 

	The expected structure of each line in the input file is as follows:
	CHROM POS ID REF ALT QUAL FILTER INFO FORMAT global ASW CEU CHB CHS CLM FIN GBR	IBS	JPT	LWK	MXL	PUR	TSI	YRI
	with each of the attributes being tab-delimited.

	The information lines explaining the data file begin with '#'; these lines are skipped when reading the file.

	The relevant aspects of each line are:
	CHROM: The chromosome number of the SNP
	POS: The position of the SNP on the chromosome
	ID: The unique ID number given to the SNP
	REF: The dominant allele for the SNP
	ALT: The recessive/alternate allele for the SNP
	global: The global expression counts for the alternate and dominant allele counts, expressed as the ratio AN:AC:SN:SC.
		AN: Total Allele Count
		AC: Alternate Allele Count
		SN: Total Sample Count*
		SC: Sample Count for each Allele*
	ASW - YRI: Expression counts for the alternate and dominant allele counts for each of the 14 individual populations on 1000 Genomes.

	*SN and SC are counts for the specific patients checked on the 1000 Genomes browser. For the purpose of these experiments, we are only concerned with AN/AC.
"""
def readPop():
	# check that the input file exists and is accessed correctly
	if os.path.isfile('populationVCF.popvcf'):
		inFile = open('populationVCF.popvcf', 'r')
	else:
		print "Could not open population file. Exiting.."
		return

	# get all lines from file
	allLines = inFile.readlines()

	# for each line in the file...
	for i in range(len(allLines)):
		# if the line begins with a '#', split to get the labels (see comments for the labels)
		if allLines[i][0] == '#':
			labels = allLines[i].split()
			# get the number of labels. this is used as a check to make sure each line of the file is not corrupted or missing information.
			numLabels = len(labels)
		else:
			# split the line by tab
			lineData = allLines[i].split('\t')
			# confirm all relevant information is in the line
			if len(lineData) == numLabels:
				# create the SNP
				makeSNPClass(lineData)

"""
	readPatient() opens the patient VCF file and reads the necessary information to create & store the individual SNP objects for each patient.
	The patient VCF file is obtained by selecting "Download Genotype Data" -> "Data: Samples" -> "Filter by: Checked Samples"

	The expected structure of each line in the input file is as follows:
	CHROM POS ID REF ALT QUAL FILTER INFO FORMAT NA19701 NA19703 NA06986 NA06989	NA18525	NA18530	HG00403	HG00404	HG01112	HG01113	HG00171	HG00182	HG00096	HG00097	HG01515	HG01516	NA18940	NA18941	NA19028	NA19036	NA19651	NA19652	HG00637	HG00638	NA20504	NA20509	NA18487	NA18489
	with each of the attributes being tab-delimited.

	The relevant aspects of each line are:
	CHROM: The chromosome number of the SNP
	POS: The position of the SNP on the chromosome
	ID: The unique ID number given to the SNP
	REF: The dominant allele for the SNP
	ALT: The recessive/alternate allele for the SNP
	NA19701 - NA18489: The 28 sample patients we checked for the larger sample.
		GT: Genotype for the patient
		DS: Genotype dosage
		GL: Genotype likelihoods

		For the purposes of these experiments, we are only concerned with the GL field.
"""
def readPatient():
	global SNPDict
	global Patients

	# check that the input file exists and is accessed correctly
	if os.path.isfile('patientVCF.vcf'):
		inFile = open('patientVCF.vcf', 'r')
	else:
		print "Could not open patient file. Exiting.."
		return

	allLines = inFile.readlines()
	labels = []

	# for each line,
	for i in range(len(allLines)):

		# get the relevant labels from the last # line.
		if allLines[i][0] == '#':
			if allLines[i+1][0] != '#':
				labels = allLines[i]
				labels = labels.strip('\n')

				# the first 9 items on each line contain information that is relevant for all patients.
				data = labels.split('\t')[0:9]

				"""
					In this method, we tested the algorithm on both a large patient file (2 samples checked for each population, for all 14 populations) and
					and a smaller patient/population combination of 3 populations, 2 samples checked for each population. Either definition of patients may
					be used, as long as it is consistent with the E assignment in initID3().
				"""
				#patients = labels.split('\t')[9:]
				patients = ['NA18487', 'NA18489', 'NA18940', 'NA18941', 'NA06986', 'NA06989']

				# Create a patient object for each patient using their ID, then add to the global list.
				for j in range(len(patients)):
					newPatient = Patient(patients[j])
					Patients.append(newPatient)

		# each line contains info about a single SNP for each of the patients
		else:
			lineData = allLines[i].split('\t')
			# As with 'patients', either if statement could be used as long as it is consistent with the rest of the experiment.
			if len(lineData) == 37:
			#if len(lineData) == 9 + len(Patients):

			# for each patient,
				for p in range(len(Patients)):
					# create a SNP using the relevant information from each line
					Patients[p].addSNP(lineData[2], lineData[p + 9], lineData[3], lineData[4])
					
					# append the patients' genotype frequencies for each SNP to a list of patientFreqs
					((SNPDict[lineData[2]]).patientFreqs).append(Patients[p].snps[lineData[2]])

				# make the average patient genotype frequency for each SNP
				SNPDict[lineData[2]].makeAvg()

"""
	Accounts for potential differences in the SNPDict and Patient SNP information - due to corruptions in either file, some patient SNPs may not be in
	the population SNPDict and vice versa.
	This method compares Patients and SNPDict, removing any SNPs that are in one but not the other.
"""
def compareSets():
	global SNPDict
	global Patients

	# get all SNP IDs from the SNPDict keys, and all SNP IDS from the keys of the first patient's SNP dictionary.
	popSNPs = SNPDict.keys()
	patSNPs = (Patients[0]).snps.keys()

	# delete SNPS from SNPDict if not also in Patients
	for p in popSNPs:
		if p not in patSNPs:
			del SNPDict[p]

	# Delete SNPs from if not in SNPDict
	for p in patSNPs:
		if p not in popSNPs:
			del SNPDict[p]	

"""
	Uses the data from each line in the populationVCF input file to construct a SNP object.
	Adds this SNP to the global dictionary.
"""
def makeSNPClass(lineData):
	global SNPDict
	
	# handles raw data file issues where lines are only partially complete by ignoring the line
	if len(lineData) == 16:
		return

	# See SNP constructor in SNPClass.py to see how each of these inputs are used.
	newSNP =  SNP(lineData[0], lineData[1], lineData[2], lineData[3], lineData[4], lineData[9], \
		lineData[10], lineData[11], lineData[12], lineData[13], lineData[14],\
		lineData[15], lineData[16], lineData[17], lineData[18], lineData[19],\
		lineData[20], lineData[21], lineData[22], lineData[23])

	# adds SNP to SNPDict using key = SNP ID, value = SNP object
	SNPDict[lineData[2]] = newSNP

"""
	Handles initialziation for the recursive ID3 algorithm. 

	In this method, E is the possible populations a patient could be classified into.
	The resulting tree will have 1 leaf per population. When the tree is traversed, the leaf node at which each patient ends up is their population
	assigned by the tree.

	F is all possible SNPs we could use to partition the different countries. The root and each internal node of the tree will be a SNP, and populations
	will be divided based on their expression ratios of the dominant and allele frequencies.

	initID3 calls ID3Recursive, which returns the root and populates the global list of all nodes (including leaves) in the decision tree.
	These are used together in traverseRec(), which traverses the decision tree and prints the classification of each sample patient tested.
"""
def initID3():
	global SNPDict 
	global Patients
	global nodes

	""" As with 'patients', either E definition could be used as long as it is consistent with the patient group being tested.
		For example, testing on the smaller 6 patient group means that E should be defined as the 3 possible populations these
		patients belong to. """
	E = ['CEU', 'JPT', 'YRI']
	#E = ['AWS', 'CEU', 'CHB', 'CHS', 'CLM', 'FIN', 'GBR', 'IBS', 'JPT', 'LWK', 'MXL', 'PUR', 'TSI', 'YRI']
	F = SNPDict.keys()

	# call the recursive ID3 algorithm
	root = ID3Recursive(SNPDict, E, F, nodes)

	# for each sample patient, traverse the tree. add the result to a dictionary storing all results in the form of:
	# key = patientID, value = classified population
	patientDict = {}
	for p in Patients:
		country = traverseRec(p, root)
		patientDict[p.id] = country

	# print the results to the user
	print  "\n\tPatient Classification Results: \n"
	for i in patientDict.items():
		print "Patient ID: ", i[0], "\tPredicted Population: ", i[1]
	print

# Main function; responsible for calling all helper methods.
def main():
	readPop()
	readPatient()
	compareSets()
	initID3()

main()