import sys
sys.path.append('../../Funcs/')
from uncModel import uncModel
import cobra
import numpy
from matplotlib import pyplot as plt
import matplotlib

# Cobrapy model to use for simulations
model = cobra.io.read_sbml_model('../../GEMs/ecYeast8-0.365.xml')
print 'Model read successfully'

# simulation of minimal media on these carbon sources, assuming glucose is set to be the only (unlimited) carbon source.
cSources = ['D-glucose exchange (reversible)', 'sucrose exchange (reversible)', 'maltose exchange (reversible)', 'D-galactose exchange (reversible)', 'ethanol exchange (reversible)', 'acetate exchange (reversible)']

# experimental values
expVal = [0.41, 0.38, 0.4, 0.28, 0.12, 0.17]  # from Van Djiken 2000, yeast strain: CEN.PK122
#exoVal = [0.44, 0.42, 0.41, 0.32, 0.12, 0.2]  # from Van Djiken 2000, yeast strain: CBS8066

# colour
colours = ['darkblue', 'darkgreen', 'lime', 'saddlebrown', 'r', 'orange']

# settings for solving the models
con = 20
Su = 1 * 10 ** -2
settingsCer = {'pathToFile': '../../GEMs/',
                'solveOutput': False,
                'Crossover': 3,
                'maxPercUncK': False,
                'linearKUncertainty': False,
                'maxKUnc': 0.2,
                'minKUnc': 0.1,
                'maxPercUncS': False,  # Su,
                'M': False,  # con * Su,
                'LPmethod': 2,
                'QPtol': 10 ** -4,
                'BarrierCorrections': -1,
                'retrieveInstance': False,
                'print': False,
                'Aggregate': 1}

settingsUnc = {'pathToFile': '../../GEMs/',
                'solveOutput': False,
                'Crossover': 3,
                'maxPercUncK': False,
                'linearKUncertainty': True,
                'maxKUnc': 0.2,
                'minKUnc': 0.1,
                'maxPercUncS': Su,   #1.51*10**-5,  # 4.75*10**-3,  # 5*10**-6 # 5*10**-4
                'M': con*Su,   #0.001,  # 1*10**-1,  # 10**-4 # 0.012
                'LPmethod': 2,
                'QPtol': 10 ** -4,
                'BarrierCorrections': -1,
                'retrieveInstance': False,
                'print': False,
                'Aggregate': 1}

# make list of settings to iterate over
settingsList = [settingsCer, settingsUnc]
settingsTitles = ['Certain model', 'Uncertain model']

# find the flux indices of interest and arrange experimental values accordingly
print('Looking for fluxes')
fluxInd = []
fluxNam = []
expValOrd = []
coloursOrd = []
allIn = 0
for i in range(len(model.reactions)):
    for j in range(len(cSources)):
        if cSources[j] == model.reactions[i].name:
            print '\tfound "' + cSources[j] + '" as "' + str(model.reactions[i].name) + '" at index ' + str(i)
            fluxInd.append(i)
            fluxNam.append(cSources[j])
            expValOrd.append(expVal[j])
            coloursOrd.append(colours[j])
            allIn += 1
if allIn == len(cSources):
    print '\tAll carbon source fluxes found'
else:
    print '\tNot all carbon source fluxes found'
print '\n',

# set all carbon source uptakes to zero
for i in range(len(fluxInd)):
    model.reactions[fluxInd[i]].upper_bound = 0

# set pyplot
fig, ax = plt.subplots(1, len(settingsList))

# iterate through different carbon sources and optimize model
for s in range(len(settingsList)):
    # print progress to screen
    print 'Starting simulations with setting ' + str(s+1) + ' of ' + str(len(settingsList)) + ' in total'

    # initialization
    grwthRates = []
    for i in range(len(fluxInd)):
        # print progress to screen
        print '\tSolving model limited by ' + fluxNam[i]

        # set given carbon source to inf
        model.reactions[fluxInd[i]].upper_bound = numpy.Inf

        # solve model
        sol, stat = uncModel(model, settingsList[s])
        if stat:
            print 'Error occured in uncModel at carbon source ' + fluxNam[i]
        grwthRates.append(sol['grwthRate'])

        # set given carbon source back to zero
        model.reactions[fluxInd[i]].upper_bound = 0
    # plot results
    ax[s].title.set_text(settingsTitles[s])
    for k in range(len(grwthRates)):
        ax[s].scatter(expValOrd[k], grwthRates[k], c=coloursOrd[k], label=fluxNam[k])
    ax[s].plot([0, 0.5], [0, 0.5])

    # print aesthetically
    print '\n',

# plot results
ax[0].legend(loc='upper left')
plt.show()

