import sys
sys.path.append('../../Funcs/')
sys.path.append('../../Funcs/GeneKnockout/')
from uncModel import uncModel
from createUncModel import createUncModel
from solveUncModel import solveUncModel
from findAllGeneKOs import findAllGeneKOs
from singleKOList import singleKOList
import cobra
import numpy
import pickle

# Simulates gene knockouts for different standard deviations of coefficients in the biomass pseudoreaction.
# Each coefficients is modeled as the only uncertain coefficient.

# Initializations of parameters
coeffStart = 0.03
coeffStop = 0.03
coeffStep = 10**1
errTol = 0.5

# initialize growth rate result list
grwthRates = []

# cobrapy model to use for simulations
model = cobra.io.read_sbml_model('../../GEMs/ecYeast8_batch.xml')
print 'Model read successfully. Model is of type:'+str(type(model))

# settings
settings = {'print': False}

# initialize RAMPER model
RAMPERpars, stat = createUncModel(model, settings)

# find growth rate for comparison
sol, stat = solveUncModel(RAMPERpars, settings)
if stat == 0:
	refGrwthRate = sol['grwthRate']
	print 'Reference growth rate = '+str(refGrwthRate)+' g/gDW'
	tolGrwthRate = (1-errTol)*refGrwthRate
	grwthRates.append(tolGrwthRate)
else:
	print 'Could not succesfully compute reference growth rate'

# find gene-reaction single knockout map
geneRxnList = singleKOList(model)

# find lethal genes
print 'Computing lethal genes for certain model'
orLethalGenes = findAllGeneKOs(model, settings, tolGrwthRate, geneRxnList)
#orLethalGenes = []
# flux names of interet
fluxNames=['pseudoreaction']

# find the flux indices of interest
print '\nLooking for biomass components'
fluxInd = []
fluxNam = []
countBioComp = 0
for i in range(len(model.reactions)):
    for j in range(len(fluxNames)):
        if fluxNames[j] in model.reactions[i].name:
            print '\t found "' + fluxNames[j] + '" as "' + str(model.reactions[i].name) + '" at index ' + str(i)
            fluxInd.append(i)
            fluxNam.append(fluxNames[j])
print 'Done searching biomass components\n'
print 'Commencing iteration through biomass components'

# model settings
settings = {'linearKUncertainty': False,
            'print': False,
	    'M': 0.0125}

# initialize results lists which will be enlargened in the next loop
resMap = ['no uncertainty']
resCoeff = ['resCoeff = '+str(coeffStart), 'M = '+str(settings['M'])]
results = [orLethalGenes]
# construct uncertainty coefficient matrix
S = cobra.util.array.create_stoichiometric_matrix(model)
# counter for coefficients. Used in results
i = 0

# start knockouts
for colIndex in fluxInd:
	reaction = S[:,colIndex]
        rowIndeces = numpy.nonzero(reaction)[0]
	print '\n\nIterating through reaction: '+str(model.reactions[colIndex].name)
	for rowIndex in rowIndeces:
		if not ((reaction[rowIndex] % 1) == 0):  # checks if coefficient is non-integer
			# update i
			i += 1
						
			# save index to resMap
			resMap.append([rowIndex, colIndex])

			# reset uncertainty matrix
			uncMat = numpy.zeros((len(model.metabolites),len(model.reactions))) 
			coeff = coeffStart

			# append results (hard-coded)
			results.append([ [], [], [] ])
			resultsInd = 0
			
			# update to screen
			print 'Metabolite '+str(model.metabolites[rowIndex].name)
			
			while coeff <= coeffStop:
				# set uncertainty coefficient
				uncMat[rowIndex, colIndex]  = coeff
				settings['uncertaintyMatrix'] = uncMat
				RAMPERpars, stat = createUncModel(model, settings)
				try:
					sol, stat = solveUncModel(RAMPERpars, settings) 
					if stat == 0:
						grwthRates.append(sol['grwthRate'])
					else:
						grwthRates.append(-1)
				except:
					grwthRates.append(-1)
				# find gene KOs				
				lethalGenes = findAllGeneKOs(model, settings, tolGrwthRate, geneRxnList)
				
				# add lethal genes to results
				results[i][resultsInd] = lethalGenes
				resultsInd += 1
				print str(sol['grwthRate'])

				# update coefficient	
				coeff = coeff * coeffStep				
# save results with pickle
with open('resultsS03.pkl', 'wb') as f:
	pickle.dump(results, f)

with open('resMapS03.pkl', 'wb') as g:
	pickle.dump(resMap, g)

with open('resCoeffS03.pkl', 'wb') as h:
	pickle.dump(resCoeff, h)

with open('grwthRatesS03.pkl', 'wb') as h:
	pickle.dump(grwthRates, h)

print 'We are done.'
