"""
    Jessica Jowdy and Brooke Kelsey 
    CS68 Final Project, Fall 2014 

    A program that creates decision trees using the training set,
    traverses the tree using test set, and prunes the tree using
    the tuning set. The program also provides a statisitcal 
    analysis of the resulting tree
"""

from SNPClass import *
from math import *


def ID3Recursive2(SNPDict2, E, F, nodes, country, parent):
    """
    Recursively creates the decision tree from the training set data
    Input: 
      SNPDict2 - dictionary of all of the SNPs common in each population 
                 data set and their dominant allele values 
      E - the elements that have attributes specified by F 
          (initially the training set) 
      F - the remaining features (initially the entire gene list) 
      nodes - the list of nodes in the current tree
      country - the specific population we are making the decision tree for
      parent - the Parent node 
    Output:
      leaf/bestGene/None - a leaf if certain conditions are met, None if other 
                       conditions are met, or bestGene, a node representing the
                       best feature to split on in the remaining feature list
    """

    if len(E) == 0:
        # no examples left to split on; do nothing 
        return None

    if len(F) == 0:
        # no features left ot split on 
        # create a leaf with the majority label at it 
        total = 0; 
        for i in range(len(E)): 
            if E[i].country == country: 
                total = total + 1
        if total > len(E)/2: 
            label = country
        else:
            label = "NO" 
        leaf = Leaf(label)
        # set the parent of the leaf
        leaf.parent = parent
	nodes.append(leaf)
	return leaf

    # determine how the examples split on the bestGene 
    branches = makeBranches2(SNPDict2, E, F, nodes,country)

    if branches == None: 
        # no good remaining features to split on 
        # create a leaf with the majority label at it 
        total = 0 
	for i in range(len(E)): 
            if E[i].country == country: 
                total = total + 1
        if total > len(E)/2: 
            label = country
        else:
            label = "NO"
        leaf = Leaf(label)
        # set the parent of the leaf 
        leaf.parent = parent
	nodes.append(leaf)
	return leaf
    
    # how the examples split based on the bestGene 
    right = branches[0]
    left = branches[1]
    bestGene = branches[2]

    bestGene.parent = parent
    # recursively call ID3Recursive2 to determine the remaining nodes in the 
    # decision tree (with newly defined example lists and feature lists)
    bestGene.rightC = ID3Recursive2(SNPDict2, right, F, nodes, country,bestGene)
    bestGene.leftC = ID3Recursive2(SNPDict2, left, F, nodes, country,bestGene)
    
    # add new node to the list 
    nodes.append(bestGene)

    return bestGene


##########################################################################
def makeBranches2(SNPDict2, E, F, nodes,country):
    """
    Segregates the examples based on allele values for the best SNP 
    to split on. 
    Input: 
      SNPDict2 - dictionary of all of the SNPs common in each population 
                 data set and their dominant allele values 
      E - the current examples list 
      F - the current features list 
      nodes - the list of nodes in the current tree
      country - the specific population we are making the decision tree for
    Output: 
      [right, left, bestGene] - the examples moving to the left child and the
            right child based on their alleles for the best SNP, and the node
            bestGene which holds the current best SNP to split on
    """
    # obtains best SNP to split on using infoGain equation
    bestGene = infoGain2(SNPDict2, E, F, country)
    
    if bestGene == None: 
        # no good feature to split on
        return None 
    
    # remove entry from SNPDict
    F.remove(bestGene.snp)
    
    # set the elements that are at bestGene 
    bestGene.setElements(E)

    right = []
    left = []
    # get dominant allele from bestGene 
    dom = bestGene.dom
    dom_allele = dom + dom  
    
    for i in range(len(E)):
        if E[i].snps[bestGene.snp] == dom_allele:
            # if dominant/dominant, split left 
            left.append(E[i]) 
        else: 
            # otherwise split right 
            right.append(E[i])
    return [right, left, bestGene] 


#######################################################################
def infoGain2(SNPDict2, E, F, country):
    """
    Calculates the infoGain for each SNP and determines which SNP maximizes 
    the calculation. A node is created with this SNP. 
    Input: 
      SNPDict2 - dictionary of all of the SNPs common in each population 
                 data set and their dominant allele values 
      E - the current examples list 
      F - the current features list 
      country - the specific population we are making the decision tree for
    Ouput: 
      geneNode - the node containing the SNP that maximizes the infoGain
    """
    total = len(E)

    max_info_gain = 0.0 
    max_gene = ""
    
    for i in range(len(F)): 
        # for each SNP in the feature list 

        # get dominant allele 
        dom = SNPDict2[F[i]]
        dom_allele = dom + dom  
        # initialize counts
        totalCD = 0
        totalCR = 0
        totalND = 0
        totalNR = 0
        country_count = 0

        for j in range(len(E)):
            # determine counts: 
            #   (CD) - how many are dominant and in population 
            #   (CR) - how many are not dominant and in population
            #   (ND) - how many are dominant and not in population 
            #   (NR) - how many are not dominant and not in population
            if E[j].snps[F[i]] == dom_allele: 
                if E[j].country == country:
                    totalCD = totalCD + 1
                    country_count = country_count + 1
                else: 
                    totalND = totalND + 1
            else: 
                if E[j].country == country: 
                    totalCR = totalCR  + 1
                    country_count = country_count + 1
                else: 
                    totalNR = totalNR + 1

        totalC = country_count  # total in country
        totalN = total - country_count  # total not in country
        totalD = totalCD + totalND  # total dominant  
        totalR = totalCR + totalNR  # total recessive 

        p_country = float(totalC)/total  # prob of being in country
        np_country = float(totalN)/total # prob of not being in country
        
        # Error checking entropy to make sure there is not 
        # a math domain error
        if p_country == 0: 
            entropy = -np_country*log(np_country,2)
        elif np_country == 0: 
            entropy = -p_country*log(p_country,2)
        else:
            entropy = -p_country*log(p_country,2)-np_country*log(np_country,2)

       
        if entropy != 0: 
            
            p_dom = totalD/total  # prob of being dominant 
            p_rec = totalR/total  # prob of being recessive 

            # Error checking probability counts to make sure there is not
            # a math domain error 
            if totalD != 0: 
                py_dom = float(totalCD)/totalD
                npy_dom = float(totalND)/totalD
            else: 
                py_dom = 0.0
                npy_dom = 0.0

            if totalR != 0: 
                py_rec = float(totalCR)/totalR
                npy_rec = float((totalR-totalCR))/totalR
            else: 
                py_rec = 0.0
                npy_rec = 0.0
      
            # Calculates first part of conditional entropy (allele is 
            # dominant) 
            if p_dom == 0: 
                cond_entropy = 0
            else:
                if py_dom == 0: 
                    cond_entropy = p_dom*(-npy_dom*log(npy_dom,2))
                elif npy_dom == 0: 
                    cond_entropy = p_dom*(-py_dom*log(py_dom,2))
                else: 
                    cond_entropy = p_dom*(-py_dom*log(py_dom,2)-npy_dom*log(npy_dom,2))
       
            # Calculates second part of conditional entropy (allele is 
            # recessive) 
            if p_rec == 0: 
                cond_entropy = cond_entropy 
            else:
                if py_rec == 0: 
                    cond_entropy = cond_entropy + p_rec*(-npy_rec*log(npy_rec,2))
                elif npy_rec == 0: 
                    cond_entropy = cond_entropy + p_rec*(-py_rec*log(py_rec,2)) 
                else: 
                    cond_entropy = cond_entropy + p_rec*(-py_rec*log(py_rec,2)-npy_rec*log(npy_rec,2))

            # subtract conditional entropy from entropy to
            # obtain the information gain 
            info_gain = entropy-cond_entropy
            
            if info_gain > max_info_gain: 
                # update the maximum infoGain and save the 
                # current SNP 
                max_info_gain = info_gain 
                max_gene = F[i]
    
    if max_gene == "":
        # no good features left to split on 
        return None

    # create a new Node with the bets SNP 
    geneNode = Node(max_gene, SNPDict2[max_gene])
        
    return geneNode


####################################################################
def traverseRec2(SNPDict2, patient, currNode):
    """
    Recursively works a patient down the tree until reaching a leaf,
    at which point the label is returned 
    Input: 
      SNPDict2 - dictionary of all of the SNPs common in each population 
                 data set and their dominant allele values 
      patient - the individual whose population is being predicted 
      currNode - the current node in the tree we are on
    Output: 
      traverseRec2() - a recursive call to the children of the current node
      currNode.country - the label of the leaf node (if curr_node 
                         is a leaf) 
    """
    global nodes

    if currNode.snp == None:
        # we reached a leaf, return the label 
        return currNode.country

    val = SNPDict2[currNode.snp]
    val_allele = val + val

    if patient.snps[currNode.snp] == val_allele:
        # if the patient is dominant/dominant, go left 
        return traverseRec2(SNPDict2, patient, currNode.leftC)
    else:
        # otherwise go right 
        return traverseRec2(SNPDict2, patient, currNode.rightC)


#####################################################################
def pruneTree(SNPDict2, nodes, root, patients, country, stats):
    """
    Determines whether a pruned tree increases the accuracy of the 
    tree, and if it does, what the pruned tree should be. 
    Input: 
      SNPDict2 - dictionary of all of the SNPs common in each population 
                 data set and their dominant allele values 
      nodes - the list of nodes in the current tree 
      root - the root of the current tree  
      patients - the tuning set used for pruning 
      country - the specific population we are making the decision tree for
      stats - the accuracy, sensitivity, and specificity of the tree
    Output: 
      [best_tree, best_stats] - the list of nodes for the best tree 
          (pruned or unpruned) and the subsequent statistics 
    """
    # initialize best tree, stats, and accuracy 
    best_tree = nodes 
    best_stats = stats
    best_accuracy = stats[0]
    # initialize node_list, curr_accuracy, curr_stats 
    node_list = [] 
    curr_accuracy = 0.0
    curr_stats = []
    for i in range(len(best_tree)):
        nodes = best_tree
        new_nodes = []
        if nodes[i].isLeaf():
            # get parent of the leaf 
            cur_node = nodes[i].parent 
            yes = 0
            no = 0
            # change parent of leaf into a leaf with the 
            # majority label 
            for k in range(len(cur_node.elements)): 
                if cur_node.elements[k].country == country: 
                    yes = yes + 1
                else: 
                    no = no + 1
            if yes > no: 
                label = country 
            else: 
                label = "NO"
            new_leaf = Leaf(label) 
            new_leaf.parent = cur_node.parent 
            
            # create a new list of nodes for pruned tree 
            for p in range(len(nodes)): 
                if nodes[p].isLeaf(): 
                    if nodes[p] != nodes[i]: 
                        new_nodes.append(nodes[p])
                else: 
                    if nodes[p] != cur_node: 
                        new_nodes.append(nodes[p])
           
            new_nodes.append(new_leaf)

            newPatientDict = {}
            
            for j in range(len(patients)): 
                # traverse tree again with tuning data set 
                new_country = traverseRec2(SNPDict2, patients[j], root) 
                newPatientDict[patients[j].id] = new_country 
            # get accuracy of new tree 
            stats = getAccuracy(newPatientDict, patients, country) 
            new_accuracy = stats[0]
            if new_accuracy > curr_accuracy: 
                # replace current accuracy, node_list, stats, with new ones
                curr_accuracy = new_accuracy
                node_list = new_nodes 
                curr_stats = stats 
    if curr_accuracy > best_accuracy: 
        # replace the unpruned tree with the pruned tree 
        best_tree = node_list 
        best_accuracy = new_accuracy 
        best_stats = curr_stats
        return [best_tree, best_stats] 
    
    return None


###################################################################
def getAccuracy(PatientDict, patients, country):
    """
    Calculates the accuracy, sensitivity, and specificity of the 
    given tree by determining the number of true positives, true
    negatives, false positives, and false negatives exist 
    Input: 
      PatientDict - a dictionary with the patient ids as keys and 
                    their predicted population as the values
      patients - the list of patients to check 
      country - the specific population we are making the decision tree for
    Output: 
      [accuracy, sensitivity, specificity] - 
    """
    TP = 0 
    TN = 0
    FP = 0
    FN = 0
    for i in range(len(patients)):
        if PatientDict[patients[i].id] == country:
            if PatientDict[patients[i].id] == patients[i].country: 
                # increment true positive 
                TP = TP + 1
            else: 
                # increment false positive 
                FP = FP + 1
        if PatientDict[patients[i].id] == "NO": 
            if patients[i].country != country: 
                # increment true negative 
                TN = TN + 1 
            else:
                # increment false negative 
                FN = FN + 1

    # complete statistical analaysis 
    accuracy = float(TP + TN)/len(patients) 
    sensitivity = float(TP)/float(TP + FN)
    specificity = float(TN)/float(TN + FP)

    return [accuracy, sensitivity, specificity]

