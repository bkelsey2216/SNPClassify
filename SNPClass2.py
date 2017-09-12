"""
    Jessica Jowdy and Brooke Kelsey 
    A program containing necessary class definitions for the 
    classification algorithm
"""

from numpy import *

class Patient2():
    """
    A class for defining the patient 
    """
    def __init__(self, ID, country):
        self.id = ID
        self.snps = {}
        self.country = country

    def addSNP(self, snpID, genotype):
        """
        Adds a SNP to the patient dictionary, where 
        the key is the SNP id and the value is the 
        patient's genotype for that SNP
        Input: 
          self 
          snpID - the id of the SNP to add 
          genotype - the patient's allele value
        Output: 
          None 
        """
        self.snps[snpID] = genotype


class Node():
    """
    A class for defining a node 
    """
    def __init__(self, SNP, dom):
        self.snp = SNP
        self.dom = dom  # dominant allele 
        self.leftC = None  # left child 
        self.rightC = None # right child 
        self.elements = []  # examples at the node 
        self.parent = None
        
    def setElements(self, lst):
        """
        Sets the elements occurring at the given node 
        Input: 
          self 
          lst - a list of examples occuring at the ndoe 
        Output: 
          None 
        """
        self.elements = lst 

    def isLeaf(self):
        """
        Determines if the current node is a leaf or not 
        Input: 
          None 
        Output: 
          False 
        """
        return False


class Leaf():
    """
    A class for defining a leaf 
    """
    def __init__(self, label):
        self.snp = None
        self.country = label  # population label
        self.parent = None

    def isLeaf(self):
        """
        Determines if the current node is a leaf or not 
        Input: 
          None 
        Output: 
          True  
        """
        return True
