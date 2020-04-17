import sys
sys.path.append('../../Funcs/')
sys.path.append('../FVA/')
from uncModel import uncModel
from pyomo.environ import Objective, Constraint, minimize, maximize, value
from pyomo.opt import SolverFactory
from numpy import linspace
import matplotlib
matplotlib.use('agg')
from matplotlib import pyplot as plt
from scipy.io import loadmat
import cobra

# TODO: make sure every upperBound is set to max

# set matplotlib rc-parameters (default parameters)
matplotlib.rc('xtick', labelsize=6)
matplotlib.rc('ytick', labelsize=6)

# Cobrapy model to use for simulations
model = cobra.io.read_sbml_model('../../GEMs/ecYeast8_batch.xml')
print 'Model read successfully. Model is of type:'+str(type(model))

# growth rates to be fixed and slack to growth rate constraint
epsPerc = 0.001
nXVal = 15  # number of growth rates to be calculated
fluxNames = ['D-glucose exchange (reversible)', 'oxygen exchange (reversible)', 'ethanol exchange', 'acetate exchange', 'carbon dioxide exchange']  # assume glucose flux is indexed first

# find the flux indices of interest
print 'Looking for fluxes'
fluxInd = []
fluxNam = []
allIn = 0
for i in range(len(model.reactions)):
    for j in range(len(fluxNames)):
        if fluxNames[j] == model.reactions[i].name:
            print '\t found "' + fluxNames[j] + '" as "' + str(model.reactions[i].name) + '" at index ' + str(i)
            fluxInd.append(i)
            fluxNam.append(fluxNames[j])
            allIn += 1
if allIn == len(fluxNames):
    print 'All input fluxes found in COBRApy model'
else:
    print 'Error: Not all input fluxes were found'
print '\n',

# substrate to minimize (glucose transport has index 93 in ecYEast17, index 1841 in ecYeast_batch, index 82 in ecYeastData) Cannot find crabtree effect in ecYeast17!!
substrUp = fluxInd[fluxNam.index(fluxNames[0])]  # assuming the first flux is glucose

#con = 22.5
	#n = 1
#Su = 5 * 10 ** -1
# set different settings for different models
settingsCer = {'pathToFile': '../../GEMs/',
               'solveOutput': False,
               'Crossover': 3,
               'maxPercUncK': False,
               'linearKUncertainty': False,
               'maxKUnc': 0.1,
               'minKUnc': 0.01,
               'maxPercUncS': False,
               'M': False,
               'LPmethod': 1,
               'QPtol': 10 ** -4,
               'BarrierCorrections': -1,
               'retrieveInstance': True,
               'print': False,
               'Aggregate': 1}

con = 50
n = 1
Su = 10 * 10 ** -1
settingsKUnc = {'pathToFile': '../../GEMs/',
                'solveOutput': False,
                'Crossover': 3,
                'maxPercUncK': False,
                'linearKUncertainty': True,
                'maxKUnc': 0.2,
                'minKUnc': 0.1,
                'maxPercUncS': False,
                'M': False,
                'LPmethod': 2,
                'QPtol': 10 ** -4,
                'BarrierCorrections': -1,
                'retrieveInstance': True,
                'print': False,
                'Aggregate': 1}

con = 22.5
#n = 1
Su = 1 * 10 ** -0
settingsS1Unc = {'pathToFile': '../../GEMs/',
                'solveOutput': False,
                'Crossover': 3,
                'maxPercUncK': False,
                'linearKUncertainty': True,
                'maxKUnc': 0.2,
                'minKUnc': 0.1,
                'maxPercUncS': Su,  # (0.5**4)*0.003,
                'M': Su*con,  # (0.5**4)*0.0625,
                'LPmethod': 2,
                'QPtol': 10 ** -4,
                'BarrierCorrections': -1,
                'retrieveInstance': True,
                'print': False,
                'Aggregate': 1}
con = 22.5
#n = 1
Su = 1 * 10 ** -1
settingsS2Unc = {'pathToFile': '../../GEMs/',
                'solveOutput': False,
                'Crossover': 3,
                'maxPercUncK': False,
                'linearKUncertainty': True,
                'maxKUnc': 0.2,
                'minKUnc': 0.1,
                'maxPercUncS': Su,  # (0.5**4)*0.003,
                'M': Su*con,  # (0.5**4)*0.0625,
                'LPmethod': 2,
                'QPtol': 10 ** -4,
                'BarrierCorrections': -1,
                'retrieveInstance': True,
                'print': False,
                'Aggregate': 1}
con = 22.5
#n = 1
Su = 1 * 10 ** -2
settingsS3Unc = {'pathToFile': '../../GEMs/',
                'solveOutput': False,
                'Crossover': 3,
                'maxPercUncK': False,
                'linearKUncertainty': True,
                'maxKUnc': 0.2,
                'minKUnc': 0.1,
                'maxPercUncS': Su,  # (0.5**4)*0.003,
                'M': Su*con,  # (0.5**4)*0.0625,
                'LPmethod': 2,
                'QPtol': 10 ** -4,
                'BarrierCorrections': -1,
                'retrieveInstance': True,
                'print': False,
                'Aggregate': 1}
con = 22.5
#n = 1
Su = 1 * 10 ** -3
settingsS4Unc = {'pathToFile': '../../GEMs/',
                'solveOutput': False,
                'Crossover': 3,
                'maxPercUncK': False,
                'linearKUncertainty': True,
                'maxKUnc': 0.2,
                'minKUnc': 0.1,
                'maxPercUncS': Su,  # (0.5**4)*0.003,
                'M': Su*con,  # (0.5**4)*0.0625,
                'LPmethod': 2,
                'QPtol': 10 ** -4,
                'BarrierCorrections': -1,
                'retrieveInstance': True,
                'print': False,
                'Aggregate': 1}
con = 22.5
#n = 1
Su = 1 * 10 ** -4
settingsS5Unc = {'pathToFile': '../../GEMs/',
                'solveOutput': False,
                'Crossover': 3,
                'maxPercUncK': False,
                'linearKUncertainty': True,
                'maxKUnc': 0.2,
                'minKUnc': 0.1,
                'maxPercUncS': Su,  # (0.5**4)*0.003,
                'M': Su*con,  # (0.5**4)*0.0625,
                'LPmethod': 2,
                'QPtol': 10 ** -4,
                'BarrierCorrections': -1,
                'retrieveInstance': True,
                'print': False,
                'Aggregate': 1}

con = 22.5
#n = 1
Su = 1 * 10 ** -5
settingsS6Unc = {'pathToFile': '../../GEMs/',
                'solveOutput': False,
                'Crossover': 3,
                'maxPercUncK': False,
                'linearKUncertainty': True,
                'maxKUnc': 0.2,
                'minKUnc': 0.1,
                'maxPercUncS': Su,  # (0.5**4)*0.003,
                'M': Su*con,  # (0.5**4)*0.0625,
                'LPmethod': 2,
                'QPtol': 10 ** -4,
                'BarrierCorrections': -1,
                'retrieveInstance': True,
                'print': False,
                'Aggregate': 1}
con = 22.5
#n = 1
Su = 1 * 10 ** -6
settingsS7Unc = {'pathToFile': '../../GEMs/',
                'solveOutput': False,
                'Crossover': 3,
                'maxPercUncK': False,
                'linearKUncertainty': True,
                'maxKUnc': 0.2,
                'minKUnc': 0.1,
                'maxPercUncS': Su,  # (0.5**4)*0.003,
                'M': Su*con,  # (0.5**4)*0.0625,
                'LPmethod': 2,
                'QPtol': 10 ** -4,
                'BarrierCorrections': -1,
                'retrieveInstance': True,
                'print': False,
                'Aggregate': 1}
con = 22.5
#n = 1
Su = 1 * 10 ** -7
settingsS8Unc = {'pathToFile': '../../GEMs/',
                'solveOutput': False,
                'Crossover': 3,
                'maxPercUncK': False,
                'linearKUncertainty': True,
                'maxKUnc': 0.2,
                'minKUnc': 0.1,
                'maxPercUncS': Su,  # (0.5**4)*0.003,
                'M': Su*con,  # (0.5**4)*0.0625,
                'LPmethod': 2,
                'QPtol': 10 ** -4,
                'BarrierCorrections': -1,
                'retrieveInstance': True,
                'print': False,
                'Aggregate': 1}
con = 22.5
#n = 1
Su = 1 * 10 ** -8
settingsS9Unc = {'pathToFile': '../../GEMs/',
                'solveOutput': False,
                'Crossover': 3,
                'maxPercUncK': False,
                'linearKUncertainty': True,
                'maxKUnc': 0.2,
                'minKUnc': 0.1,
                'maxPercUncS': Su,  # (0.5**4)*0.003,
                'M': Su*con,  # (0.5**4)*0.0625,
                'LPmethod': 2,
                'QPtol': 10 ** -4,
                'BarrierCorrections': -1,
                'retrieveInstance': True,
                'print': False,
                'Aggregate': 1}
con = 22.5
#n = 1
Su = 1 * 10 ** -9
settingsS10Unc = {'pathToFile': '../../GEMs/',
                'solveOutput': False,
                'Crossover': 3,
                'maxPercUncK': False,
                'linearKUncertainty': True,
                'maxKUnc': 0.2,
                'minKUnc': 0.1,
                'maxPercUncS': Su,  # (0.5**4)*0.003,
                'M': Su*con,  # (0.5**4)*0.0625,
                'LPmethod': 2,
                'QPtol': 10 ** -4,
                'BarrierCorrections': -1,
                'retrieveInstance': True,
                'print': False,
                'Aggregate': 1}
# different models to plot against each other
modelSettings = []
#modelSettings.append(settingsCer)
#modelSettings.append(settingsKUnc)
modelSettings.append(settingsS1Unc)
modelSettings.append(settingsS2Unc)
modelSettings.append(settingsS3Unc)
modelSettings.append(settingsS4Unc)
modelSettings.append(settingsS5Unc)
#modelSettings.append(settingsS6Unc)
#modelSettings.append(settingsS7Unc)
#modelSettings.append(settingsS8Unc)

# experimental data
expGrwth = [0.025, 0.05, 0.1, 0.15, 0.2, 0.25, 0.28, 0.3, 0.35, 0.4]
expData = [[], [], [], [], []]
expData[fluxNam.index(fluxNames[0])] = [0.3, 0.6, 1.1, 1.7, 2.3, 2.8, 3.4, 4.5, 8.6, 11.1]  # glucose expenditure
expData[fluxNam.index(fluxNames[1])] = [0.8, 1.3, 2.5, 3.9, 5.3, 7, 7.4, 6.1, 5.1, 3.7]  # oxygen expenditure
expData[fluxNam.index(fluxNames[2])] = [0, 0, 0, 0, 0, 0, 0.11, 2.3, 9.5, 13.9]  # ethanol production
expData[fluxNam.index(fluxNames[3])] = [0, 0, 0, 0, 0, 0, 0.08, 0.41, 0.62, 0.6]  # acetate production
expData[fluxNam.index(fluxNames[4])] = [0.8, 1.4, 2.7, 4.2, 5.7, 7.5, 8, 8.8, 14.9, 18.9]  # carbon dioxide production

# colours and label for each model type (to be used in plotting)
#col = ['b', 'r', 'y']
#col = ['b']
#for i in range(len(modelSettings)):
#	c = next(col[i])
#	col.append(c)
#col = range(len(modelSettings))

# create labels to be used as subtitles for each model in final plot
lab = []
for i in range(len(modelSettings)):
	modelSetting = modelSettings[i]
	if modelSetting['linearKUncertainty']:
		kStr = 'K: ' + str(int(round(100*modelSetting['minKUnc']))) + '-' + str(int(round(100*modelSetting['maxKUnc']))) + '%, '
	else:
		kStr = 'K: 0%, '
	if modelSetting['maxPercUncS']:
		sStr = 'S: ' + str(100*modelSetting['maxPercUncS']) + '%\n'
	else:
		sStr = 'S: 0%\n'
	if modelSetting['M']:
		mStr = 'M: ' + str(modelSetting['M'])
	else: 
		mStr = 'M: 0'
	#subtitle = kStr + sStr + mStr
	subtitle = sStr + mStr
	lab.append(subtitle)

# settings for solving chemostat simulations
opt = SolverFactory('gurobi')
opt.options['BarHomogeneous'] = 1  # appears necessary for uncertain yeast models
opt.options['NumericFocus'] = 3  # seems to give odd optimal objectives if not set to 3
opt.options['BarQCPConvTol'] = 10**-4

# plot config.
fig, ax = plt.subplots(len(fluxInd), len(modelSettings), sharex='all', sharey='row')

# iterate over different types of models (certain, uncertain etc.)
for i in range(len(modelSettings)):
    # set settings
    settings = modelSettings[i]

    # initialize lists of fluxes to plot
    fluxPlt = [[] for _ in range(len(fluxInd))]

    # print progress to screen
    print 'Starting at model ' + str(i + 1) + ' of ' + str(len(modelSettings))

    # solve an optimization problem in order to get the concrete model
    sol, stat = uncModel(model, settings)
    instance = sol['instance']
    cS = sol['cS']
    MW = sol['MW']
    optGrwthRate = sol['grwthRate']
    instance.grwthRate.deactivate()

    # set growth rate x-values based on optimal growth rate
    grwthRates = linspace(0, optGrwthRate, nXVal).tolist()

    # initialize upper and lower bound on optimal fluxes, fva[0] - lower bounds, fva[1] - upper bounds
    fva = [[[] for _ in range(len(fluxInd))] for _ in range(2)]

    # check if growth rates are feasible in model
    if grwthRates[-1] > optGrwthRate:
        print '\tError: maximum growth rate in input exceeds optimal growth rate attainable for model: ' + str(optGrwthRate)
        exit()
    else:
        print '\tOriginal instance retrieved. Optimal growth rate: ' + str(optGrwthRate)

    # iterate over growth rates
    for j in range(len(grwthRates)):
        # set growth rate
        grwthRate = grwthRates[j]

        # alter instance to constrain biomass and minimize substrate uptake
        def bioConstrLow(instance):
            return sum(cS[j] * instance.v[j] for j in instance.flux) <= (1 + epsPerc) * grwthRate

        instance.bioConLow = Constraint(rule=bioConstrLow)
        instance.bioConLow.activate()

        def bioConstrHigh(instance):
            return sum(cS[j] * instance.v[j] for j in instance.flux) >= (1 - epsPerc) * grwthRate

        instance.bioConHigh = Constraint(rule=bioConstrHigh)
        instance.bioConHigh.activate()
        instance.substrUptake = Objective(expr=instance.v[substrUp], sense=minimize)
        instance.substrUptake.activate()

        # solve for smallest substrate uptake
        Soln = opt.solve(instance)
        minSubstr = value(instance.substrUptake)

        # deactivate objective and constaint
        instance.substrUptake.deactivate()


        # constraint substrate uptake to optimal value
        def substrConstrLow(instance):
            return instance.v[substrUp] <= (1 + epsPerc) * minSubstr


        def substrConstrHigh(instance):
            return instance.v[substrUp] >= (1 - epsPerc) * minSubstr


        instance.substrConLow = Constraint(rule=substrConstrLow)
        instance.substrConHigh = Constraint(rule=substrConstrHigh)

        # add objective minimizing total enzyme mass
        instance.TotMass = Objective(expr=sum(MW[j] * instance.e[j] for j in instance.enzK), sense=minimize)
        instance.TotMass.activate()

        # solve for smallest enzyme mass
        Soln = opt.solve(instance)

        # find flux values to plot
        for s in range(len(fluxInd)):
            fluxPlt[s].append(value(instance.v[fluxInd[s]]))

        # constrain minimal total enzyme mass and deactivate totMass objective
        minTotMass = value(instance.TotMass)
        def totMassLow(instance):
            return sum(MW[j] * instance.e[j] for j in instance.enzK) >= (1 - epsPerc) * minTotMass
        def totMassHigh(instance):
            return sum(MW[j] * instance.e[j] for j in instance.enzK) <= (1 + epsPerc) * minTotMass
        instance.totMassLow = Constraint(rule=totMassLow)
        instance.totMassHigh = Constraint(rule=totMassHigh)
        instance.TotMass.deactivate()

        # FVA analysis of fluxes
        for f in range(len(fluxInd)):
            # activate totMass constraints
            instance.totMassLow.activate()
            instance.totMassHigh.activate()

            # We set an upper bound on the flux in question to avoid unbounded problems.
            def fluxConstr(instance):
                return instance.v[fluxInd[f]] <= max(expData[f]) + 10
            instance.fluxConstr = Constraint(rule=fluxConstr)
            instance.fluxConstr.activate()

            # maximize
            instance.maxv = Objective(expr=instance.v[fluxInd[f]], sense=maximize)
            instance.maxv.activate()
            Soln = opt.solve(instance)
            maxv = value(instance.maxv)

            # minimize
            instance.maxv.deactivate()
            instance.minv = Objective(expr=instance.v[fluxInd[f]], sense=minimize)
            instance.minv.activate()
            Soln = opt.solve(instance)
            minv = value(instance.minv)
            instance.minv.deactivate()

            # store max and min values
            fva[0][f].append(fluxPlt[f][j]-minv)
            fva[1][f].append(maxv-fluxPlt[f][j])

            # deactivate
            instance.totMassHigh.deactivate()
            instance.totMassLow.deactivate()
            instance.fluxConstr.deactivate()

        # deactivate constraints and objectives
        instance.bioConLow.deactivate()
        instance.bioConHigh.deactivate()
        instance.substrConLow.deactivate()
        instance.substrConHigh.deactivate()
        instance.TotMass.deactivate()

        # print progress to screen
        print '\tSolved instance ' + str(j + 1) + ' of ' + str(len(grwthRates))

    # plot model
    # if multiple settings tested (necessary because of multi-indexing in ax does not work if only one setting is ran)
    if len(modelSettings) > 1:
        ax[0, i].set_title(lab[i], fontdict={'fontsize': '7'})
   
        ax[len(fluxInd) - 1, (len(modelSettings)//2)].set_xlabel('Growth rate (g/gDW)', fontsize=8, labelpad=5)       # 0.5125, 0.02, 'Growth rate (g/gDWh)', horizontalalignment='center', fontsize=8)
	for l in range(len(fluxInd)):		
            ax[l, i].errorbar(grwthRates, fluxPlt[l], yerr=[fva[0][l], fva[1][l]], fmt='o', markersize=2)  # yerr=[fva[0][l] for _ in range(len(fluxPlt[l]))])
            ax[l, i].scatter(expGrwth, expData[l], marker='v', c='g', s=3)

        # set y- and x-axis labels
        for l in range(len(fluxInd)):
            if fluxNam[l][-1] == ')':
                fixedFluxNam = fluxNam[l][0:-13]
            else:
                fixedFluxNam = fluxNam[l]
            ax[l, 0].set_ylabel(fixedFluxNam[0:-9] + '\n' + fixedFluxNam[-8:], fontsize=7)
    else:
        ax[0].title.set_text(lab[i])
        ax[len(fluxInd) - 1].set_xlabel('Growth rate (g/gDW)', fontsize=8, labelpad=5)       # 0.5125, 0.02, 'Growth rate (g/gDWh)', horizontalalignment='center', fontsize=8)
        for l in range(len(fluxInd)):
            ax[l].scatter(expGrwth, expData[l], marker='v', c='g')
            ax[l].errorbar(grwthRates, fluxPlt[l], yerr=[fva[0][l], fva[1][l]], fmt='o')  # yerr=[fva[0][l] for _ in range(len(fluxPlt[l]))])

        # set y- and x-axis labels
        for l in range(len(fluxInd)):
            if fluxNam[l][-1] == ')':
                fixedFluxNam = fluxNam[l][0:-13]
            else:
                fixedFluxNam = fluxNam[l]
            ax[l].set_ylabel(fixedFluxNam[0:-9] + '\n' + fixedFluxNam[-8:])

# show/save figure (has to check for number of settings to index ax correctly)
#if len(modelSettings) > 1:
#    ax[0, 0].legend(loc='upper left')
#else:
#    ax[0].legend(loc='upper left')
plt.savefig('./Results/crabecYeast8_batch_K10-20.pdf')  # subplots overlap, you may avoid this by saving the matlab figure as a pdf instead
#plt.show()

# we are done
print '\nDone.'
