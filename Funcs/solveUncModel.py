from uncUtils import *
import numpy
import time

# TODO: Make way for an option to only retrieve the pyomo istance without solving the model
def solveUncModel(parDict, inpSettings={}):
    # assume failure
    sol = declareNone(1)
    stat = 1
    
    # flag for suboptimal solution
    subopt = False
    
    # Set default settings if none is chosen
    settings = {'logFile': 'logFile',
                'seeWarnings': False,
                'solver': 'gurobi',
                'LPmethod': 2,
                'NumericFocus': 3,
                'QPtol': 10 ** (-8),
                'LPtol': 10 ** (-8),
                'Crossover': 0,
                'BarrierCorrections': -1,
                'solveOutput': False,
                'pathToFile': '',
                'maxPercUncK': False,
                'maxPercUncS': False,
                'M': False,
                'matModName': '',
                'print': True,
                'retrieveInstance': False,
                'linearKUncertainty': False,
                'maxKUnc': 0.2,
                'minKUnc': 0.05,
                'percKForced': [1],
                'Presolve': 2,
                'BarHomogeneous': -1,
		'uncertaintyMatrix': [],
		'covarianceList': [],
                'Aggregate': 1}  # will be properly initialized later.


    # remove old log file
    if os.path.exists(settings['logFile']):
        os.remove(settings['logFile'])

    # Change settings according to input
    if set(inpSettings.keys()).issubset(settings.keys()):
        for key, item in inpSettings.items():
            settings[key] = item
    else:
        mesg(settings['logFile'], 'Error: could not read variables in createUncModel()\n')
        print 'Error: check logFile'
        return sol, stat

    # check if keys are correct in input
    keysColl = {'SsprCer': 0,
		  'SsprUnc': 0,
	          'Kspr': 0,
		  'P': 0,
 		  'MW': 0,
		  'bSCer': 0,
		  'bSUnc': 0,
		  'bK': 0,
		  'cS': 0,
  	          'ubS': 0,
 		  'lbS': 0,
		  'ubK': 0,
 		  'lbK': 0,
  		  'totMass': 0,
 		  'RSspr': 0,
		  'RSindx': 0,
 		  'RKspr': 0,
	          'RKindx': 0,
 		  'M': 0, 
		  'cerRowInd': 0,
 		  'uncerRowInd': 0,
 		  'stat': 0}
    
    
    if not set(parDict.keys()).issubset(keysColl.keys()):
    	mesg(settings['logFile'], 'Input of solveUncModel is incomplete with respect to parameters.\n') 	

    try:
        # Solve model
        sol, solveStat = solveModel(parDict['SsprCer'], parDict['SsprUnc'], parDict['Kspr'], parDict['P'], parDict['MW'], parDict['bSCer'], parDict['bSUnc'], parDict['bK'], parDict['cS'], parDict['ubS'], parDict['lbS'], parDict['ubK'], parDict['lbK'], parDict['totMass'], parDict['RSspr'], parDict['RSindx'], parDict['RKspr'], parDict['RKindx'], parDict['M'], parDict['cerRowInd'], parDict['uncerRowInd'],settings)
        if solveStat == 1:
            print 'Error occurred in solveModel(). Read logFile for more information'
            return sol, solveStat
        if solveStat == 2:
            subopt = True
        mesg(settings['logFile'], 'Evaluated modelSolve() successfully in uncModel()\n')
    except:
        mesg(settings['logFile'], 'Error: could not evaluate modelSolve()\n')
        print 'Error: check logFile'
        return sol, stat
    # we are done, print growth rate if told to
    stat = 0
    if settings['print']:
        if subopt:
            print 'Model solved suboptimally.\n\nSuboptimal growth rate = '+str(sol['grwthRate'])
        else:
            print 'Model solved successfully.\n\nOptimal growth rate = '+str(sol['grwthRate'])
   
    # final print to screen
    if settings['print']:
	print 'Done.\n' 

    return sol, stat
