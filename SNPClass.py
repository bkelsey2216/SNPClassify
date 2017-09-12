from numpy import *

"""
	A class to store all information relevant to a single SNP.
	Has 3 methods: 
		makeAvg() - make the average genotype frequencies using all sample patient frequencies
		getCountry() - return the country object, given the country label.
		getAllCountries() - return a tuple of all country objects
"""
class SNP():

	""" all constructor parameters are read in from the population VCF file. this currently creates attributes for every possible
	# population in the 1000 Genomes browser. to consider fewer populations per SNP, this will have to be modified to be compatible with
	# any modifications to the input file format. """
	def __init__(self, chro, pos, ID, refAll, altAll, Globalinfo, AWSinfo, CEUinfo, CHBinfo, CHSinfo, CLMinfo,\
	FINinfo, GBRinfo, IBSinfo, JPTinfo, LWKinfo, MXLinfo, PURinfo, TSIinfo, YRIinfo):
		# position on the gene
		self.position = pos
		# chromosome number
		self.chromosome = chro
		# unique ID of the SNP
		self.id = ID
		# country objects for the global counts and 14 individual countries
		self.Global = Country(Globalinfo)
		self.AWS = Country(AWSinfo)
		self.CEU = Country(CEUinfo)
		self.CHB = Country(CHBinfo)
		self.CHS = Country(CHSinfo)
		self.CLM = Country(CLMinfo)
		self.FIN = Country(FINinfo)
		self.GBR = Country(GBRinfo)
		self.IBS = Country(IBSinfo)
		self.JPT = Country(JPTinfo)
		self.LWK = Country(LWKinfo)
		self.MXL = Country(MXLinfo)
		self.PUR = Country(PURinfo)
		self.TSI = Country(TSIinfo)
		self.YRI = Country(YRIinfo)
		# alternate and dominant (ref) alleles
		self.refAll = refAll
		self.altAll = altAll
		# list of patient genotype frequencies
		self.patientFreqs = []
		# average of all  patient genotype frequencies
		self.freqMed = [0.0, 0.0, 0.0]

	# once patientFreqs is full, construct the freqMed tuple by taking
	# the median of patientFreqs
	def makeAvg(self):
		total = [[], [], []]
		for i in range(len(self.patientFreqs)):
			for j in range(0, 3):
				total[j].append(self.patientFreqs[i][j])

		for i in range(0,3):
			self.freqMed[i] = median(total[i])

	# return country object based on string label passed in
	def getCountry(self, label):
		if label == "AWS":
			return self.AWS
		if label == "CEU":
			return self.CEU 
		if label == "CHB":
			return self.CHB 
		if label == "CHS":
			return self.CHS 
		if label == "CLM":
			return self.CLM 
		if label == "FIN":
			return self.FIN 
		if label == "GBR":
			return self.GBR 
		if label == "IBS":
			return self.IBS 
		if label == "JPT":
			return self.JPT 
		if label == "LWK":
			return self.LWK 
		if label == "MXL":
			return self.MXL 
		if label == "PUR":
			return self.PUR 
		if label == "TSI":
			return self.TSI 
		if label == "YRI":
			return self.YRI

	# return tuple of all country objects
	def getAllCountries(self):
		countries = [self.AWS, self.CEU, self.CHB, self.CHS, self.CLM, self.FIN, self.GBR,\
		 self.IBS, self.JPT, self.LWK, self.MXL, self.PUR, self.TSI, self.YRI]
		return countries

"""
	A class to store the allele counts for each population, taken from population VCF. These values are referenced
	when calculated information gain, and are stored as attributes of the SNP class.
"""
class Country():

	def __init__(self, info):
		infoSplit = info.split(":")
		# total allele count
		self.AN = float(infoSplit[0])
		# alternate allele count
		self.AC = float(infoSplit[1])
		# sample count for each allele
		self.SN = float(infoSplit[2])
		# total sample count
		self.SC = float(infoSplit[3])

"""
	A patient class, storing all SNP data for an individual patient and their patient ID.
	Method addSNP() handles converting the genotype frequencies (dominant-dominant, dominant-recessive, recessive-recessive)
	from log probabilities into base-10 probabilities that sum to 1. It then adds the SNP information to the patient's SNP dictionary
	with key = SNPID, value = 3 tuple of genotype frequencies.
"""
class Patient():
	def __init__(self, ID):
		self.id = ID
		self.snps = {}

	def addSNP(self, snpID, genotype, Dom, Rec):
		genotype = genotype.strip()
		genoSplit = genotype.split(':')
		freqs = genoSplit[2].split(',')

		AA = 10**(float(freqs[0]))
		Aa = 10**(float(freqs[1]))
		aa = 10**(float(freqs[2]))

		self.snps[snpID] = [AA, Aa, aa]

"""
	A node class used to store the SNP information of an internal node in the decision tree. Also stores
	the threshold used to partition populations, and the left and right children (leaves or nodes) as
	computed by ID3Recursive.
"""
class Node():

	def __init__(self, SNP):
		self.snp = SNP
		self.threshold = 0.0
		self.leftC = None
		self.rightC = None

"""
	A leaf class storing the population label for each SNP, and also a snp field set to None
	which is used to differentiate between internal nodes and leaves as we traverse the tree to classify a patient.
"""
class Leaf():

	def __init__(self, label):
		self.snp = None
		self.country = label