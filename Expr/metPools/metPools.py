from Funcs.uncModel import uncModel
from scipy.io import loadmat
import numpy as np
from matplotlib import pyplot as plt
import warnings
import matplotlib.cbook
warnings.filterwarnings("ignore", category=matplotlib.cbook.mplDeprecation)

# path to data
dataPath = 'C:/Users/cga32/OneDrive/Masteroppgave/PythonScriptsandFiles/RAMPGECKO/GEMs/'

# values to plot
uncSList = [10 ** -4, 10 ** 0, 10 ** 4]
Mlist = [50, 55, 60, 65, 70]
metsToPlot = np.array(range(92))
stopName = 10  # stops naming metabolites

# Set input to solver
model = 'ecYeast_batchUb.mat'
settings = {'pathToFile': dataPath,
            'Crossover': 3,
            'maxPercUncK': False,
            'maxPercUncS': 0.1,
            'M': 1000,
            'LPmethod': 2,
            'solveOutput': False,
            'print': False}

# Find names of metabolites
pathToFile = settings['pathToFile']
matFile = pathToFile + model
modInfo = loadmat(matFile)
matModName = [modInfo.keys()[i] for i in range(len(modInfo.keys())) if modInfo.keys()[i] != '__version__' and
              modInfo.keys()[i] != '__header__' and modInfo.keys()[i] != '__globals__'][0]
metName = modInfo[matModName]['metNames']  # list of metabolite ( + protein) names
metName = np.array([str(metName[0][0][j][0][0]) for j in range(len(metName[0][0]))])

# Create plots of change in concentration over different uncertain parameters
nOfProb = len(Mlist)*len(uncSList)
fig, ax = plt.subplots(1, len(Mlist), sharex='col', sharey='row')
for j in range(len(Mlist)):
    settings['M'] = Mlist[j]
    for i in range(0, len(uncSList)):
        settings['maxPercUncS'] = uncSList[i]
        # solve model, look in uncUtils for structure of sol.
        sol, stat = uncModel(model, settings)
        print 'Solved problem '+str(i + 1 + len(uncSList)*j)+' of '+str(nOfProb)
        if stat:
            print 'Error: could not execute uncModel(). See logFile'
        else:
            # Find change in concentration of metabolites
            v = np.array(sol['optFlux'])
            SsprUnc = sol['SsprUnc']
            metPools = SsprUnc.dot(v)
            # Find metabolites that we want to plot
            metPools = metPools[metsToPlot]
            uncerRowInd = np.array(sol['uncerRowIndS'])
            metNameij = metName[uncerRowInd][metsToPlot]
            #metName = metName[metsToPlot]
            # plot
            if j == 0:
                ax[j].set_ylabel(r'Expected change in concentration')
                if len(metsToPlot) < stopName:
                    plt.setp(ax, xticks=range(len(metPools)), xticklabels=metName)
            if len(metsToPlot) < stopName:
                #ax[j].set_xticks(range(len(metPools)), metName)
                #plt.setp(ax, xticks=range(len(metPools)), xticklabels=metName)
                plt.sca(ax[j])
                plt.xticks(range(len(metPools)), metName, rotation='vertical', fontsize=5)
            else:
                ax[j].set_xticks([])
            #plt.sca(ax[j])
            #plt.xticks(range(len(metPools)), metName, rotation='vertical', fontsize=5)
            ax[j].scatter(range(len(metPools)), metPools, s=4)
            ax[j].set_title(r'M = '+str(Mlist[j]))
            #plt.axhline(y=0, linestyle='--')  # horizontal line through zero


#plt.show()
fig.legend(['Perc. Unc. in S = ' + str(uncSList[i]) for i in range(len(uncSList))], loc='upper center', ncol=len(Mlist))
fig.savefig('metPools2.pdf', bbox_inches='tight')






