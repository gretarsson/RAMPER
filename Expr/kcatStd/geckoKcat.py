import sys
sys.path.append('../../Funcs/')
from uncModel import uncModel
from createUncModel import createUncModel
from solveUncModel import solveUncModel
import cobra
import numpy as np
import pickle
import time
from math import exp
from functools import partial
from multiprocessing import Pool

# Here we increase each individual standard deviation to infer the sensitivity of increasing an enzyme's turnover numbers.
# We also model each kcat per the regression model extrapolated from BRENDA data.

# function for incrementally increasing enzyme's turnover number. Formulated as to enable multiprocessing
def findMaxStd(stoichMat, grwthCutoff, tnCons, maxIter, protRows, i):
	from math import exp

	# initialize
	results = [[], [], []]	

	# load model anew 
	model = cobra.io.read_sbml_model('../../GEMs/ecYeast8_batch.xml')
	
	# Get the K matrix
	settings = {'print': False}
	RAMPERpars, stat = createUncModel(model, settings)
	Kspr = RAMPERpars['Kspr']

	# find enzyme row and its turnover number column indices
	rowInd = i
	
	# shallowly copy the original matrix row that is to be configured
	orRow = Kspr[rowInd]	

	# change each turnover number (non-integer values) simultaneously, and save indeces with kcats for later
	kcatVar = []
	for j in range(orRow.shape[1]):
		if not orRow[0,j].is_integer():
			Kspr[rowInd,j] = tnCons[0]*orRow[0,j]
			kcatVar.append(j)

	# save info on this particular comp. experiment.
	results[0].append(stoichMat[protRows[rowInd]])
	results[1].append(protRows[rowInd])
	
	# calculate growth rate
	counter = 1
	try:
		settings = {'print': False}
		RAMPERpars['Kspr'] = Kspr
		RAMPERsol, stat = solveUncModel(RAMPERpars, settings)
	except:
		stat = 1
	if stat:	
		results[2].append('error')
		tnGrwthRate = grwthCutoff + 1  # so that we do not stop iterating
	else:	
		tnGrwthRate = RAMPERsol['grwthRate']
		
		# save results
		results[2].append(tnGrwthRate)

	# update counter and reset the matrix row (this is to be certain of the multiplier)
	counter += 1	

	while tnGrwthRate > grwthCutoff and counter <= maxIter:
		Kspr[rowInd,kcatVar] = tnCons[counter-1]*orRow[0,kcatVar]
		try:
			settings = {'print': False}
			RAMPERpars['Kspr'] = Kspr
			RAMPERsol, stat = solveUncModel(RAMPERpars, settings)
		except:
			stat = 1
		if stat:	
			results[2].append('error')
			tnGrwthRate = grwthCutoff + 1  # so that we do not stop iterating
		else:
			tnGrwthRate = RAMPERsol['grwthRate']
				
			# save results
			results[2].append(tnGrwthRate)

		# update counter and reset row
		counter += 1

	# print to screen
	print '\tCompleted simulations for enzyme '+str(i+1)+' of total '+str(len(protRows))
	return results	


# ------------------------------------------------------------------------------------------------------
def main():
	from math import exp
	
	# time script
	start_time = time.time()	

	# parameters for comp. expr.
	grwthCutoff = 0.5  # portion for which we say model is unstable
	maxIter = 10  # max iteration of increase in enzyme's kcats before we stop
	#tnCons = np.linspace(0, 0.01, maxIter)
	tnCons = [10**(i+1) for i in range(maxIter)]

	# initializations of results
	results = []  # final standard deviation as the whole, final stoichmetric uncertainty matrix row,  amount of iterations. final growth rate, corresponding row index in GEM for above results

	# use regression model to create stoichiometric uncertainty matrix
	modelFile = '../../GEMs/ecYeast8_batch.xml'
	model = cobra.io.read_sbml_model(modelFile)
	stoichMat = cobra.util.array.create_stoichiometric_matrix(model)
	print  'Model read successfully'

	protRows = [] # save row indices for enzyme constraints
	# iterate metabolites and find protein rows
	for j in range(len(model.metabolites)-1): # avoid last index i.e. the pool constraint
		rowName = model.metabolites[j].name
		if 'prot_' == rowName[0:5]:
			# save row index
			protRows.append(j)		

	# Iterate through individual turnover numbers and increase their standard deviation.
	settings = {'print': False}
	try:
		RAMPERpars, stat = createUncModel(model, settings)
		RAMPERsol, stat = solveUncModel(RAMPERpars, settings)
	except:
		print 'Could not compute reference growth rate'
	grwthRate = RAMPERsol['grwthRate']
	grwthCutoff = grwthCutoff * grwthRate
	print 'Computed reference growth rate successfully: '+str(grwthRate)+' g/gDWh-1'

	# start multiprocessing
	pool = Pool(18)
	func = partial(findMaxStd, stoichMat, grwthCutoff, tnCons, maxIter, protRows) 
	allResults = pool.map(func, range(len(protRows))) 
	print 'Completed all simulations successfully'
	metaData = {	'tnCons': tnCons,
			'grwthCutoff': grwthCutoff,
			'maxIter': maxIter,
			'Model name': modelFile,
			'note': 'No uncertainty, tnCons increases the inversed turnover number (decreases turnover) of a GECKO version of the model'
		    } 
	# dump the data using pickle
	outAll = open('kcatExpGECKO.p', 'wb')
	outMeta = open('metaKcatExpGECKO.p', 'wb')
	pickle.dump(allResults, outAll)
	pickle.dump(metaData, outMeta)
	outAll.close()
	outMeta.close()
	print 'Saved pickle files successfully'
	print '\nScript executed in '+str(time.time() - start_time)+' seconds.'

	# we are done.

if __name__ == '__main__':
	main()


