import sys
sys.path.append('../')
from uncModel import uncModel
from createUncModel import createUncModel
from solveUncModel import solveUncModel
import cobra
import numpy

# Inputs a cobra model (model) with its settings (settings) and iteratively finds all lethal gene knockouts below a user-set growth rate cutoff (grwthCut) solved using the uncModel function
# The output consists of a list the size of genes. The elements of the list has three values: -1 (could not solve model), 0 (growth rate not lethal), and 1 (lethal).
def findAllGeneKOs(model, settings, grwthCut, geneRxnList, compareList=[]):
	# create uncertainty model
	RAMPERpars, stat = createUncModel(model, settings)
	lbS = RAMPERpars['lbS']
	ubS = RAMPERpars['ubS']
	lbK = RAMPERpars['lbK']
	ubK = RAMPERpars['ubK']
	colSepInd = len(lbS)	
	
	# solve model for reference growth rate
	try:
		sol, stat = solveUncModel(RAMPERpars, settings)
		grwthCut = sol['grwthRate']	
	except:
		pass

        # extract number of genes
	nOfGenes = len(geneRxnList)
	lethalGenes = [-1 for ele in range(nOfGenes)]
	
	# instantiate index map for cobrapy reactions in rxnGeneList
	r_ind = model.reactions.index	

	# iterate through the reactions (genes)
	for i in range(nOfGenes):
		# identify reactions to be blocked and initialize
		reactions = geneRxnList[i]
		lb = [0 for reaction in reactions]
		ub = [0 for reaction in reactions]
		r_indices = [0 for reaction in reactions]		

		# change and save default reaction bounds
		for j in range(len(reactions)):
			# find index in stoichimetric matrix
			r_indices[j] = r_ind(reactions[j])
			stoichInd = r_indices[j]
			# Checks whether reaction is in S or K 
			if stoichInd < colSepInd:
				# save default bounds
				lb[j] = lbS[stoichInd]
				ub[j] = ubS[stoichInd]
				
				# knockout reaction/gene
				RAMPERpars['lbS'][stoichInd]  = 0
				RAMPERpars['ubS'][stoichInd]  = 0
			else:
				# save default bounds
				lb[j] = lbK[stoichInd-colSepInd]
				ub[j] = ubK[stoichInd-colSepInd]
				
				# knockout reaction/gene
				RAMPERpars['lbK'][stoichInd-colSepInd]  = 0
				RAMPERpars['ubK'][stoichInd-colSepInd]  = 0
				
		# optimize using uncModel
		try:
			sol, stat = solveUncModel(RAMPERpars, settings)
			if not stat:
				grwthRate = sol['grwthRate']
				lethalGenes[i] = int(grwthRate < grwthCut)
				if compareList:
					if compareList[i] != lethalGenes[i]:
						lethalGenes.append('Ended abruptly.')
		except:
			pass
	
		# set bounds back to normal
		for j in range(len(reactions)):
			stoichInd = r_indices[j]
			if stoichInd < colSepInd:
				RAMPERpars['lbS'][stoichInd] = lb[j]
				RAMPERpars['ubS'][stoichInd] = ub[j]
			else:
				RAMPERpars['lbK'][stoichInd-colSepInd] = lb[j]
				RAMPERpars['ubK'][stoichInd-colSepInd] = ub[j]
						
	# we are done
	return lethalGenes
	
