import cobra
import sys
sys.path.append('../../Funcs/')
from uncModel import uncModel
from createUncModel import createUncModel
from solveUncModel import solveUncModel
import time
import numpy as np
from math import exp

# In this script, we compute the average computation time of RAMPER (with regression model) and ecFBA problems as performed by the RAMPER implementation.
# We will test the computation time of reading models/parsing data, constructing Pyomo model instances, and solving the optimization problems.

# initializations
n_list = [100]  # number of times each computation time is performed
model = cobra.io.read_sbml_model('../../GEMs/ecYeast8_batch.xml')  # read COBRApy model

# Compute deviation matrix of the regression model
# constants in regression model
cons = exp(-0.73)
m = 1.1

# extract stoichiometric matrix from COBRApy model and instantiate deviation matrix
stoichMat = cobra.util.array.create_stoichiometric_matrix(model)
uncMat = np.zeros((len(model.metabolites), len(model.reactions)))

# iterate metabolites and find protein rows
for j in range(len(model.metabolites) - 1):  # avoid last index i.e. the pool constraint
    rowName = model.metabolites[j].name
    if 'prot_' == rowName[0:5]:
        # save row index
        row = stoichMat[j, :]
        colInds = np.nonzero(row)[0]
        # iterate through turnover numbers and find standard deviation via regression model
        for colInd in colInds:
            if row[colInd].is_integer():
                continue  # skip iteration if coefficient is not a turnover number
            std = cons * (abs(row[colInd]) ** m)
            # put standard deviation into stoichiometric uncertainty matrix
            uncMat[j, colInd] = std

# start computations
for n in n_list:
    # instantiate lists to save results
    ecFBA_createUnc_time = []  # time for ecFBA in createUncModel()
    ecFBA_solveUnc_time = []  # time for ecFBA in solveUncModel()
    ecFBA_solver_time = []  # time for Gurobi to solve Pyomo model instance
    RAMPER_createUnc_time = []  # analogous for RAMPER
    RAMPER_solveUnc_time = []  # analogous for RAMPER
    RAMPER_solver_time = []  # analogous for RAMPER

    #print 'Setting n = '+str(n)
    # ecFBA
    settings = {'print': False}
    #print 'Running ecFBA...'
    for i in range(n):
        start = time.time()
        RAMPERpars, stat = createUncModel(model, settings)
        end = time.time()
        ecFBA_createUnc_time.append(end-start)
        del stat, start, end

        start = time.time()
        sol, stat = solveUncModel(RAMPERpars, settings)
        end = time.time()
        ecFBA_solveUnc_time.append(end-start)
        ecFBA_solver_time.append(sol['solverTime'])
        del sol, stat, RAMPERpars, start, end

        #print '\tCompleted '+str(i+1)+' of '+str(n)

    # RAMPER

    # put deviation matrix into settings
    settings = {'uncertaintyMatrix': uncMat,
                'print': False}

    #print '\nRunning RAMPER...'
    for i in range(n):
        start = time.time()
        RAMPERpars, stat = createUncModel(model, settings)
        end = time.time()
        RAMPER_createUnc_time.append(end - start)
        del stat, start, end

        start = time.time()
        sol, stat = solveUncModel(RAMPERpars, settings)
        end = time.time()
        RAMPER_solveUnc_time.append(end - start)
        RAMPER_solver_time.append(sol['solverTime'])
        del sol, stat, RAMPERpars, start, end

        #print '\tCompleted ' + str(i+1) + ' of ' + str(n)

    # print results to screen
    print '----------------------------------'
    print 'Results for n = '+str(n)
    print 'ecFBA computation time averages:'
    print '\tcreateUnc: '+str(sum(ecFBA_createUnc_time)/len(ecFBA_createUnc_time))
    print '\tsolveUnc: '+str(sum(ecFBA_solveUnc_time)/len(ecFBA_solveUnc_time))
    print '\tInstiating Pyomo: '+str(sum(ecFBA_solveUnc_time)/len(ecFBA_solveUnc_time) - sum(ecFBA_solver_time)/len(ecFBA_solver_time))
    print '\tGurobi solver: '+str(sum(ecFBA_solver_time)/len(ecFBA_solver_time))


    print '\nRAMPER computation time averages:'
    print '\tcreateUnc: '+str(sum(RAMPER_createUnc_time)/len(RAMPER_createUnc_time))
    print '\tsolveUnc: '+str(sum(RAMPER_solveUnc_time)/len(RAMPER_solveUnc_time))
    print '\tInstiating Pyomo: '+str(sum(ecFBA_solveUnc_time)/len(ecFBA_solveUnc_time) - sum(ecFBA_solver_time)/len(ecFBA_solver_time))
    print '\tGurobi solver: '+str(sum(RAMPER_solver_time)/len(RAMPER_solver_time))
    print '----------------------------------'

