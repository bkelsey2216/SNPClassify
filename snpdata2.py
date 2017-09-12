"""
    Jessica Jowdy and Brooke Kelsey 
    CS68 Final Project, Fall 2014 

    A program that runs classification on a set of patient data 
    to determime which population a patient belongs to. The 
    program also validates the results and provides a statistical
    analysis of its accuracy. 
"""


from SNPClass2 import *
from DTree2 import *
from random import *
import os 

SNPDict = {}
Patients = []
test = []
SNPDict2 = {}

def main():
        gene_list = []
        genes = []
        total_patients = []
        files = ['/local/YRI_2.txt', '/local/CEU_2.txt', '/local/JPT_2.txt']
        country = ["YRI", "CEU", "JPT"]
        gene_list = getSNPs(files[0],0,gene_list)
        for i in range(len(files)-1): 
            gene_list = getSNPs(files[i+1],1,gene_list)
        for j in range(len(files)): 
            total_patients.append(readPatient2(files[j], country[j], gene_list))
            
        data = crossValidation(total_patients)
	initID3_2(data, country, gene_list)


#################################################################3#
def getSNPs(filename, flag, old_genes):
    """
    A program that compares the SNPs in each of the patient files
    and creates a list of only the ones common in each of the 
    population data files 
    Input: 
      filename - the name of the patient file 
      flag - determines whether we have an initial gene list or not
      old_genes - the gene list from the previous comparison of 
                  patient files
    Output:
      old_genes - the common genes between the first i data sets in 
                  the filename list
    """
    new_genes = []
    dele = []
    inFile = open(filename, 'r')
    allLines = inFile.readlines()
    for i in range(len(allLines)):
        if i != 0:
            # add the SNP id to the population (new) gene list
            lineData = allLines[i].split(' ')
            new_genes.append(lineData[0])
   
    if flag == 0: 
        # create initial gene list of the population SNPs
        old_genes = new_genes
            
    else: 
        # compare population gene list to existing gene list 
        for k in range(len(old_genes)):
            if old_genes[k] not in new_genes: 
                # if element in existing gene list is not 
                # in the new population gene list, add it to 
                # list of elements to remove 
                dele.append(old_genes[k])

        for j in range(len(dele)):
            # remove each element from the existing gene list 
            old_genes.remove(dele[j])

    # return the revised gene list 
    return old_genes


######################################################################
def readPatient2(filename, country, genes):
    """
    Creates a patient class for each of the patients in the 
    population data set and saves the allele value of the 
    patient for that particular SNP in the patient's SNP 
    dictionary. The program also creates a global dictionary
    that stores the SNP's dominant allele. 
    Input: 
      filename - the name of the patient file 
      country - the population for the patient file 
      genes - the list of common SNPs among each population data
              file 
    Output: 
      country_patients - the list of patients in the population 
                         data set 
    """
    country_patients = []
    
    if os.path.isfile(filename):
       inFile = open(filename, 'r')
    else:
        print "File not valid" 
        return 
    
    allLines = inFile.readlines()
    for i in range(len(allLines)):
        if i == 0:
            # create a Patient with the patient ID and the
            # population they belong to 
	    labels = allLines[i]
	    labels = labels.strip('\n')
	    data = labels.split(' ')[0:11]
	    patients = labels.split(' ')[11:]
	    for j in range(len(patients)):
                newPatient = Patient2(patients[j],country)
	        country_patients.append(newPatient)
                    
	else:
	    lineData = allLines[i].split(' ')
            if lineData[0] in genes: 
                # if the current SNP is in the gene list 
                # add the patient's allele value for that
                # particular SNP into its personal SNP Dictionary
	        for p in range(len(patients)):
	            country_patients[p].addSNP(lineData[0], lineData[p + 11])
                # add dominant allele for SNP to SNP Dictionary 
                # for easy access later 
                SNPDict2[lineData[0]] = lineData[1][0]
        
    return country_patients


#####################################################################
def crossValidation(data):
    """ 
    Creates the training and test sets for the program. The
    training set will be used to create the decision tree and 
    the test set will be ued to evaluate the accuracy of the tree
    Input: 
      data - a list of patient population lists 
    Output: 
      [training_set, test_set, tuning_set] - the training, test, and 
                         tuning sets to be used on the ID3 algorithm
    """
    training_set = []
    test_set = []
    tuning_set = []
    
    for i in range(len(data)): 
        while len(data[i]) > 10:
            # randomly choose elements in each patient population 
            # list to be in the training set 
            x = randrange(len(data[0]))
            training_set.append(data[i][x])
            # remove the elements from the list 
            del data[i][x]

        while len(data[i]) > 5:
            # randomly choose elements in each patient population
            # list to be in tuning set 
            x = randrange(len(data[0]))
            tuning_set.append(data[i][x]) 
            # remove the elements from the list
            del data[i][x]

    for j in range(len(data)):
        for k in range(len(data[j])):
            # add the five remaining patients from each population
            # into the test set 
            test_set.append(data[j][k])

    return [training_set, test_set, tuning_set] 

#####################################################################
def initID3_2(data, country_list, gene_list):
    """
    A program that initializes the ID3 Algorithm, creates the decision
    tree model, and validates the model with the test set data. The 
    program also runs prunning on the existing tree until a (local) 
    maximum accuracy is reached 
    Input:
      data - the training and test sets  
      country_list - the list of populations 
      gene_list - the list of SNPs common in the populations 
    Output: 
      None
    """
    global SNPDict2
    # training set is example list in ID3 Algorithm 
    E = data[0]
    # SNP data is feature list in ID3 Algorithmn 
    F = gene_list
    test = data[1]
    tuning = data[2]
    nodes = []
    # make the decision tree using the training data  
    root = ID3Recursive2(SNPDict2, E, F, nodes, country_list[2],None)
    
    patientDict = {}
    for i in range(len(test)):
        # validate the tree using the test set data 
        country = traverseRec2(SNPDict2, test[i], root)
        patientDict[test[i].id] = country
        print test[i].id, test[i].country, country
    
    # get initial accuracy using test set data 
    stats = getAccuracy(patientDict, test, country_list[2])
    while pruneTree(SNPDict2, nodes, root, tuning, country_list[2], stats) != None: 
        # prune tree until we have reach a (local) maximum accuracy 
        info = pruneTree(SNPDict2, nodes, root, tuning, country_list[2], stats)
        nodes = info[0]
        stats = info[1]
    print stats





main()
