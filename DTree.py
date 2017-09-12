"""
	Jessica Jowdy, Brooke Kelsey
	CS68 Final Project, Fall 2014

	Library containing all functions for the first decision tree method.
"""

from SNPClass import *
from math import *

"""
	ID3Recursive() performs the recursive step of the ID3 algorithm, partitioning the examples (populations) down a binary tree until
	leaves have been created for every feature. The root of this tree is returned to initID3().
"""
def ID3Recursive(SNPDict, E, F, nodes):

	# if there are no populations left to partition, return None
	if len(E) == 0:
		return None

	# if there are no SNPs left to use as internal nodes or only one population was sent down the branch,
	if len(F) == 0 or len(E) == 1:
		# create a leaf object for that population,
		leaf = Leaf(E[0])
		# append it to the global list of tree nodes, and return the leaf 
		nodes.append(leaf)
		return leaf
	
	# call makeBranches(), a function that calls the infoGain() to select the best SNP to partition on
	# and partitions E into left and right branches
	branches = makeBranches(SNPDict, E, F, nodes)

	""" If the left or right branch has no populations, or there is no genetic variance between any of the populations
	# (an average of '1.0' indicates that all patients definitely have the dominant allele, so splitting at that allele
	# gives us no information), continue calling makeBranches()"""
	while len(branches[0]) == 0 or len(branches[1]) == 0 or branches[2].snp.freqMed[0] == 1.0:
		branches = makeBranches(SNPDict, E, F, nodes)

	# branches is a 3 tuple of [right populations, left populations, partitioning SNP]. extract these values:
	right = branches[0]
	left = branches[1]
	bestGene = branches[2]
	
	# append the partitioning SNP to the tree nodes list as an internal node/root node
	nodes.append(bestGene)

	# call ID3Recursive() on the left and right branches. The return values from these functions
	# will be the left and right children (both are attributes of the Node class) of the current partitioning SNP.
	bestGene.rightC = ID3Recursive(SNPDict, right, F, nodes)
	bestGene.leftC = ID3Recursive(SNPDict, left, F, nodes)

	# additionally set the left and right 'edges' as the alternate and dominant alleles. The resulting tree
	# is structured such that the populations most dominant for the allele are further left.
	bestGene.rightE = bestGene.snp.altAll
	bestGene.leftE = bestGene.snp.refAll

	# return the node up the tree. at the top level, this returns the root.
	return bestGene

"""
	makeBranches() tests each remaining SNP for the amount of information it provides to partition populations.
	A SNP where all populations have identical gene expression rates is not an informative gene. It selects the best SNP to partition
	on, then uses the population's aggregate counts to create left and right branches. It returns all of this information
	to ID3Recursive().
"""
def makeBranches(SNPDict, E, F, nodes):

	# Call infoGain to locate the most informative SNP. bestGene is of type Node (see SNPClass.py)
	bestGene = infoGain(SNPDict, F, E)
	

	# remove entry for this SNP from SNPDict, so that we don't try to use it as another internal node
	del SNPDict[bestGene.snp.id]

	# remove entry from F
	F.remove(bestGene.snp.id)

	# get the recessive allele expression ratio for each population
	ratios = []
	for i in range(len(E)):
		# get the country object (see SNPClass.py) for the particular population and particular SNP
		country = bestGene.snp.getCountry(E[i])
		# add the recessive/total ratio to the list of all ratios
		ratios.append(country.AC/country.AN)

	# set the threshold for partitioning using the median of these ratios
	avgRatio = median(ratios)
	bestGene.threshold = avgRatio

	# divide the current list of populations, E, into left and right branches based on their recessive/total ratio
	# compared to the threshold ratio for the gene.
	right = []
	left = []
	for i in range(len(E)):
		country = bestGene.snp.getCountry(E[i])
		countryRatio = country.AC/country.AN

		if countryRatio < avgRatio:
			right.append(E[i])
		else:
			left.append(E[i])

	# return a tuple containting [list of right populations, list of left populations, SNP to create internal node for]
	return [right, left, bestGene] 

"""
	infoGain() is the information gain function central to feature selection in the ID3 algorithm. In this implementation,
	it selects the best SNP by calulating the effect that partitioning based on a given SNP would have on the overall entropy
	for each population. The information gain for each individual population is calculated, and these gains are summed.
	
	The SNP providing the highest aggregate information gain is chosen as the best node to partition on.
	A Node object is created for the SNP and returned to makeBranches(). 
"""
def infoGain(SNPDict, F, pops):

	# initialize the maximum information gain to be 0.
	IGmax = [0.0, 'cat']
	
	# for each of the features,
	for i in range(len(F)):

		# totalD = total number of patients with the dominant allele
		totalD = SNPDict[F[i]].Global.AN - SNPDict[F[i]].Global.AC
		# totalR = total number of patients with the recessive allele
		totalR = SNPDict[F[i]].Global.AC
		# global aggregate = sum of these two values
		totalDR = totalD + totalR

		# initialize the entropy of the SNP using the global count information
		# if there are no patients with the dominant allele,
		if totalD == 0:
			# all terms using totalD in the numerator are ignored when calculating the initial entropy
			HY = - (totalR/totalDR)*log(totalR/totalDR, 2)
		# if there are no patients with the recessive allele,
		elif totalR == 0:
			# all terms using totalR in the numerator are ignored
			HY = -(totalD/totalDR)*log(totalD/totalDR, 2)
		# else, both are involved in the entropy calculation.
		else:
			HY = -(totalD/totalDR)*log(totalD/totalDR, 2) - (totalR/totalDR)*log(totalR/totalDR, 2)

		IG = 0.0
		# for each population,
		for j in range(len(pops)):

			# get the Country object stored in the SNP
			country = SNPDict[F[i]].getCountry(pops[j])

			# compute the population counts, just like we did for the global counts
			totalCD = country.AN - country.AC
			totalCR = country.AC
			totalC = totalCD + totalCR

			# aggregate information gain into IG

			# if no patients in the population have the dominant allele,
			if totalCD == 0:
				# calculate information gain
				IG = IG - (totalC/totalDR)*(totalCR/totalC)*log(totalCR/totalC, 2)
			# if no patients in the population have the recessive allele,
			elif totalCR == 0:
				# calculate the information gain
				IG = IG - (totalC/totalDR)*(totalCD/totalC)*log(totalCD/totalC, 2)
			# else, use both to calculate the gain
			else:
				IG = IG - (totalC/totalDR)*((totalCD/totalC)*log(totalCD/totalC, 2) + (totalCR/totalC)*log(totalCR/totalC, 2))

		# if the (intial entropy - aggregate info gain) is greater than the current maximum information gain, replace
		# the current maximum and the best SNP
		if (HY-IG) >= IGmax[0]:
		    IGmax = [IG, F[i]]

	# create node for the best SNP
	geneNode = Node(SNPDict[IGmax[1]])

	return geneNode

"""
	Recursively traverse the tree using the patient's information and the average dominant/recessive frequency for each SNP.
	Once at a leaf, return the predicted population to initID3().
"""
def traverseRec(patient, currNode):
	global nodes

	# If at a leaf, return the population label stored at the leaf.
	if currNode.snp == None:
		return currNode.country 

	# If the patient's recessive allele frequency at the current SNP node is less than the SNP's average, go left
	if currNode.snp.freqMed[0] > patient.snps[currNode.snp.id][0]:
		return traverseRec(patient, currNode.leftC)
	# otherwise, go right
	else:
		return traverseRec(patient, currNode.rightC)