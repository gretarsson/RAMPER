from uncUtils import *
import numpy


# TODO: Make way for an option to only retrieve the pyomo istance without solving the model
def uncModel(modelFile, inpSettings):
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
        mesg(settings['logFile'], 'Error: could not read variables in uncModel()\n')
        print 'Error: check logFile'
        return sol, stat

    # set bool variable for uncertainty matrix 
    uncMat = numpy.array(settings['uncertaintyMatrix'])
    uncMatBool = uncMat.any()

    # write to logFile/screen
    mesg(settings['logFile'], 'Running uncModel()\n')
    if settings['print']:
        print 'Initiating modelling'

    try:
        # Parse data
        SsprCer, Kspr, P, MW, bSCer, bK, cS, ubS, lbS, ubK, lbK, totMass, slmInd, Kmap, Smap, SmodMap, parseStat = parseData(modelFile, settings)
        if parseStat:
            print 'Error occurred in parseData(). Read logFile for more information'
            return sol, stat
        mesg(settings['logFile'], 'Evaluated parseData() successfully in uncModel()\n\n')
        if settings['print']:
            print 'Parsing data completed'
    except:
        mesg(settings['logFile'], 'Error: Could not evaluate parseData() in uncModel()\n')
        print 'Error: check logFile'
        return sol, stat

    try:
        # Model K uncertainty or not
        if bool(settings['maxPercUncK']) != bool(settings['linearKUncertainty']):
            # initialize
            if settings['linearKUncertainty'] and len(settings['percKForced']) == 1:
                settings['percKForced'] = [settings['percKForced'][0]]*Kspr.shape[0]
            # model K uncertainty
            RKspr, RKindx, mKstat = modelKcatUncertainty(Kspr, Kmap, settings)
            if mKstat:
                print 'Error occurred in modelKcatUncertainty(). Read logFile for more information'
                return sol, stat
            if settings['print']:
                print 'Modelling K uncertainty completed'
            mesg(settings['logFile'], 'Evaluated modelKcatUncertainty() successfully in uncModel()\n\n')
        elif bool(settings['maxPercUncK']) or bool(settings['linearKUncertainty']):
            print 'Error: check logFile'
            mesg(settings['logFile'], 'Error: set either maxPercUncK or linearKuncertainty to False in uncModel()\n')
            return sol, stat
        else:
            # set K as certain
            RKspr = []
            RKindx = []
            mesg(settings['logFile'], 'No K uncertainty input: Skipped modelKcatUncertainty() in uncModel()\n')
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
                return sol, stat
            # sets default value if not set by user
            if not settings['M']:
                settings['M'] = 0
            if settings['print']:
                print 'Modelling S uncertainty completed'
            mesg(settings['logFile'], 'Evaluated modelSUncertainty() successfully in uncModel()\n\n')
        else:
            # set S as certain
            SsprUnc = []
            bSUnc = []
            RSspr = []
            RSindx = []
            cerRowInd = []
            uncerRowInd = []
            settings['M'] = []
            mesg(settings['logFile'], 'No S uncertainty input: Skipped modelSUncertainty() in uncModel()\n\n')
    except:
        mesg(settings['logFile'], '\nError: could not evaluate modelSUncertainty()\n\n')
        print 'Error: check logFile'
        return sol, stat

    try:
        # Set M
        settings['M'] = [settings['M']] * SsprCer.shape[0]
    except:
        mesg(settings['logFile'], 'Error: could not calculate M vector in uncModel()\n')
        print 'Error: check logFile'
        return sol, stat

    try:
        # Solve model
        sol, solveStat = solveModel(SsprCer, SsprUnc, Kspr, P, MW, bSCer, bSUnc, bK, cS, ubS, lbS, ubK, lbK, totMass, RSspr, RSindx, RKspr, RKindx, settings['M'], cerRowInd, uncerRowInd,settings)
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
