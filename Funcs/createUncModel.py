from uncUtils import *
import numpy


# TODO: Make way for an option to only retrieve the pyomo istance without solving the model
def createUncModel(modelFile, inpSettings={}):
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
        return outputDict, stat

    # set bool variable for uncertainty matrix 
    uncMat = numpy.array(settings['uncertaintyMatrix'])
    uncMatBool = uncMat.any()

    # write to logFile/screen
    mesg(settings['logFile'], 'Running createUncModel()\n')
    if settings['print']:
        print 'Initiating modelling'

    try:
        # Parse data
        SsprCer, Kspr, P, MW, bSCer, bK, cS, ubS, lbS, ubK, lbK, totMass, slmInd, Kmap, Smap, SmodMap, parseStat = parseData(modelFile, settings)
        if parseStat:
            print 'Error occurred in parseData(). Read logFile for more information'
            return outputDict, stat
        mesg(settings['logFile'], 'Evaluated parseData() successfully in createUncModel()\n\n')
        if settings['print']:
            print 'Parsing data completed'
    except:
        mesg(settings['logFile'], 'Error: Could not evaluate parseData() in createUModel()\n')
        print 'Error: check logFile'
        return outputDict, stat

    try:
        # Model K uncertainty or not
        if bool(settings['maxPercUncK']) != bool(settings['linearKUncertainty']) or uncMatBool:
            # initialize
            if settings['linearKUncertainty'] and len(settings['percKForced']) == 1:
                settings['percKForced'] = [settings['percKForced'][0]]*Kspr.shape[0]
            # model K uncertainty
            RKspr, RKindx, mKstat = modelKcatUncertainty(Kspr, Kmap, settings)
            if mKstat:
                print 'Error occurred in modelKcatUncertainty(). Read logFile for more information'
                return outputDict, stat
            if settings['print']:
                print 'Modelling K uncertainty completed'
            mesg(settings['logFile'], 'Evaluated modelKcatUncertainty() successfully in uncModel()\n\n')
        elif bool(settings['maxPercUncK']) or bool(settings['linearKUncertainty']):
            print 'Error: check logFile'
            mesg(settings['logFile'], 'Error: set either maxPercUncK or linearKuncertainty to False in createUncModel()\n')
            return outputDict, stat
        else:
            # set K as certain
            RKspr = []
            RKindx = []
            mesg(settings['logFile'], 'No K uncertainty input: Skipped modelKcatUncertainty() in createUncModel()\n')
    except:
        print 'Error: check logFile'
        mesg(settings['logFile'], '\nError: could not evaluate modelKcatUncertainty()\n')

    try:
        # Model S unertainty or not
        if settings['maxPercUncS'] or settings['M'] or uncMatBool:
            # sets default value if not set by user
            if not settings['maxPercUncS']:
                settings['maxPerUncS'] = 0
            # model S uncertainty
            SsprCer, SsprUnc, bSCer, bSUnc, RSspr, RSindx, cerRowInd, uncerRowInd, mSstat = modelSUncertainty(SsprCer, bSCer, slmInd, Smap, settings)
            if mSstat:
                print 'Error occurred in modelSUncertainty(). Read logFile for more information'
                return outputDict, stat
            # sets default value if not set by user
            if not settings['M']:
                settings['M'] = 0
            if settings['print']:
                print 'Modelling S uncertainty completed'
            mesg(settings['logFile'], 'Evaluated modelSUncertainty() successfully in createUncModel()\n\n')
        else:
            # set S as certain
            SsprUnc = []
            bSUnc = []
            RSspr = []
            RSindx = []
            cerRowInd = []
            uncerRowInd = []
            settings['M'] = []
            mesg(settings['logFile'], 'No S uncertainty input: Skipped modelSUncertainty() in createUncModel()\n\n')
    except:
        mesg(settings['logFile'], '\nError: could not evaluate modelSUncertainty()\n\n')
        print 'Error: check logFile'
        return outputDict, stat

    try:
        # Set M
        settings['M'] = [settings['M']] * SsprCer.shape[0]
    except:
        mesg(settings['logFile'], 'Error: could not calculate M vector in createUncModel()\n')
        print 'Error: check logFile'
        return outputDict, stat

    
    stat = 0
    mesg(settings['logFile'], 'Done creating uncertainty model\n')
    outputDict = {'SsprCer': SsprCer,
		  'SsprUnc': SsprUnc,
	          'Kspr': Kspr,
		  'P': P,
 		  'MW': MW,
		  'bSCer': bSCer,
		  'bSUnc': bSUnc,
		  'bK': bK,
		  'cS': cS,
  	          'ubS': ubS,
 		  'lbS': lbS,
		  'ubK': ubK,
 		  'lbK': lbK,
  		  'totMass': totMass,
 		  'RSspr': RSspr,
		  'RSindx': RSindx,
 		  'RKspr': RKspr,
	          'RKindx': RKindx,
 		  'M': settings['M'], 
		  'cerRowInd': cerRowInd,
 		  'uncerRowInd': uncerRowInd,
		  'SmodMap': SmodMap,
 		  'stat': stat}
    return outputDict, stat
