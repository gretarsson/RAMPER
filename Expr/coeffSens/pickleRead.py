import pickle
import cobra

with open('./results/resultsS02.pkl', 'rb') as f:
	results = pickle.load(f)
with open('./results/resMapS02.pkl', 'rb') as f:
	resMap = pickle.load(f)
with open('./results/resCoeffS02.pkl', 'rb') as f:
	resCoeff = pickle.load(f)
with open('./results/grwthRatesS02.pkl', 'rb') as f:
	grwthRates = pickle.load(f)

# read cobra model
model = cobra.io.read_sbml_model('../../GEMs/ecYeast8_batch.xml')
print 'Done reading cobra model.\n'

# find resCoeff that gives different resultS
orResult = results[0]
biocoeffRes = []
# iterate through biomass coefficients
for i in range(1,len(results)):
	resultCoeff = results[i]
	index = resMap[i][0]
	metName = model.metabolites[index].name 
	grwthRate = grwthRates[i-1]
	print str(metName)+', growth rate: '+str(grwthRate)
	# iterate through different uncertainty values
	for j in range(len(resCoeff)):
		result = resultCoeff[j]
		diffIndices = [s for s in range(len(result)) if result[s] != orResult[s] and result[s] != -1 and orResult[s] != -1]
		biocoeffRes.append(resCoeff[j])
		printRes = '\t'+str(resCoeff[j]) + ': ' + str(len(diffIndices))
		print printRes
	print '\n'
print 'Done'
