from Funcs.uncModel import uncModel
from pyomo.environ import Constraint
from pyomo.opt import SolverFactory
import pickle
from multiprocess import Pool, cpu_count
from functools import partial
import time
import os
# TODO: gurobi is hardcoded as solver (true of all code in RAMPER really),
#       fix assumption that v >= 0 for all v,
#       round the growth rate in logFile to first four nonzero digits


def FVA(model, settings, fileName='', fluxes=-1, inpFvaSettings=None, epsPerc=0.01, lowerBound=10 ** -8, upperBound=10 ** 3, nOfProcessors=cpu_count() / 2):
    # Functions have to be included above parallelization task. Cannot use import in normal sense. Windows only probably.
    def mesg(logFile, mesgToAdd):
        fl = open(logFile, "a+")
        fl.write(mesgToAdd)
        fl.close()

    def fvaSubProblem(Orinstance, OrOpt, size, optV, lowerBound, upperBound, i):
        # import statements need to be inside due to parallelization, only necessary for Windows probably
        from pyomo.environ import Objective, minimize, maximize, value, Constraint
        from pyomo.opt import TerminationCondition
        import logging

        # Suppress pyomo warnings
        logging.getLogger('pyomo.core').setLevel(logging.ERROR)

        #  initialize (min, max, termination condition at min, termination condition at max, index, opt. v w/ biomass)
        fva = [[], [], [], [], [], []]

        # (ABSOLUTELY NECESSARY) though not entirely sure why...
        instance = Orinstance
        Opt = OrOpt

        # minimize flux
        if optV[i] < lowerBound:
            minState = 'Flux is lower than' + str(lowerBound) + 'in optimal solution'
            minV = 0
        else:
            try:
                # change objective to minimizing flux
                instance.minFlux = Objective(expr=instance.v[i], sense=minimize)  # default sense is minimize
                instance.minFlux.activate()

                # try to avoid sub. term.
                instance.boxMin = Constraint(expr=instance.v[i] >= lowerBound)
                instance.boxMin.activate()

                # solve min. problem
                Soln = Opt.solve(instance, tee=False)
                if Soln.solver.termination_condition == TerminationCondition.unbounded:
                    minV = lowerBound
                elif Soln.solver.termination_condition == TerminationCondition.optimal:
                    minV = value(instance.minFlux)
                else:
                    minV = -1
                minState = str(Soln.solver.termination_condition)
                instance.minFlux.deactivate()
                instance.boxMin.deactivate()
            except Exception as e:
                instance.minFlux.deactivate()
                instance.boxMin.deactivate()
                minState = 'Cannot load a SolverResults object with bad status: error'
                minV = -1
                # print(e)

        # maximize flux
        if optV[i] >= upperBound:
            maxState = 'Flux is larger than ' + str(upperBound) + ' in optimal solution.'
            maxV = upperBound
        else:
            try:
                # change objective to maximizing flux
                instance.maxFlux = Objective(expr=instance.v[i], sense=maximize)
                instance.maxFlux.activate()

                # set box constraint on flux to roughly predict unboundedness.
                instance.boxMax = Constraint(expr=instance.v[i] <= upperBound)
                instance.boxMax.activate()

                # solve max. problem
                Soln = Opt.solve(instance, tee=False)
                if Soln.solver.termination_condition == TerminationCondition.unbounded:
                    maxV = upperBound
                elif Soln.solver.termination_condition == TerminationCondition.optimal:
                    maxV = value(instance.maxFlux)
                else:
                    maxV = -1
                maxState = str(Soln.solver.termination_condition)
                instance.maxFlux.deactivate()
                instance.boxMax.deactivate()
            except Exception as e:
                instance.maxFlux.deactivate()
                instance.boxMax.deactivate()
                maxState = 'Cannot load a SolverResults object with bad status: error'
                maxV = -1
                # print(e)

        # update FVA list
        fva[0] = minV
        fva[1] = maxV
        fva[2] = minState
        fva[3] = maxState
        fva[4] = i
        fva[5] = optV[i]

        # check for errors to log
        minError = ''
        maxError = ''
        if minV == -1:
            minError = '. Error in minimization problem: ' + minState
        if maxV == -1:
            maxError = 'Error in maximization problem: ' + maxState

        # send update to screen
        mesg(settings['logFile'], '\tSolved FVA subproblem ' + str(i) + ' of ' + str(size - 1) + minError + '. ' + maxError + '\n')

        # we're done
        return fva

    # following if-sentence necessary in Windows for use of multiprocess package
    if __name__ == '__main__':
        # initialize
        fvaData = []

        # Sets logfile to default filename and checks if logFile already exists
        if 'logFile' not in settings:
            settings['logFile'] = 'logFileFVA'
        if os.path.exists(settings['logFile']):
            os.remove(settings['logFile'])

        try:
            # set fvaSettings to default if none is chosen
            fvaSettings = {'solver': 'gurobi', 'NumericFocus': 3, 'Aggregate': 0, 'BarQCPConvTol': 10 ** -3,
                           'BarHomogeneous': 1, 'Presolve': -1, 'DualReductions': 1, 'BarCorrectors': -1, 'Method': 2}
            if inpFvaSettings is not None:
                if set(inpFvaSettings.keys()).issubset(fvaSettings.keys()):
                    for key, item in inpFvaSettings.items():
                        fvaSettings[key] = item
        except:
            mesg(settings['logFile'], 'Error: could not load settings argument properly')
            return fvaData


        mesg(settings['logFile'], 'Initializing FVA process\n\n')

        # measure time
        start_time = time.time()

        # set solver options for FVAsubproblem. Interior point experiences less problems differentiating feasibility/unboundedness
        # setting all settings to default except for Aggregate, which is turned off (=0), and BarQCPConvTol set to 10**-3, seems to work fine.
        solver = fvaSettings['solver']
        OrOpt = SolverFactory(solver)
        OrOpt.options['NumericFocus'] = fvaSettings['NumericFocus']  # 0 default. When set to zero, multiple max. problems cannot differentiate infeasible and unbounded problems.
        OrOpt.options['Presolve'] = fvaSettings['Presolve']  # -1 default
        OrOpt.options['Aggregate'] = fvaSettings['Aggregate']  # Presolver options. Gurobi docs recommend turning this off if num. trouble, 1 - default, 0 - off
        OrOpt.options['DualReductions'] = fvaSettings['DualReductions']  # put to zero to verify if unbounded or infeasible
        OrOpt.options['BarQCPConvTol'] = fvaSettings['BarQCPConvTol']
        OrOpt.options['BarHomogeneous'] = fvaSettings['BarHomogeneous']  # -1 default (automatic), 0 off, 1 on: Some min problems results in unloadable SolveResults object if not turned on
        OrOpt.options['BarCorrectors'] = fvaSettings['BarCorrectors']
        OrOpt.options['Method'] = fvaSettings['Method']

        # set this automatically since we will need the original instance anyway
        settings['retrieveInstance'] = True

        mesg(settings['logFile'], 'Starting to solve initial optimization problem to find optimal growth rate\n')
        # fetch opt. growth rate and flux state and pyomo concrete model.
        try:
            sol, stat = uncModel(model, settings)
            grwthRate = sol['grwthRate']
            Oinstance = sol['instance']
            optV = sol['optFlux']
            cS = sol['cS']
        except:
            mesg(settings['logFile'], 'Error: could not solve initial optimization problem\n')
            return fvaData
        mesg(settings['logFile'], 'Initial optimization problem solved. Optimal growth rate is '+str(grwthRate)+' g/gDWh\n\n')

        # constrain biomass reaction
        def bioConstr(instance):
            return sum(cS[j] * instance.v[j] for j in instance.flux) >= (1 - epsPerc) * grwthRate

        # fix growth rate and deactivate biomass as objective
        Oinstance.con = Constraint(rule=bioConstr)
        Oinstance.grwthRate.deactivate()

        # set flux to default
        if type(fluxes) is not list and fluxes == -1:
            fluxes = range(len(optV))

        mesg(settings['logFile'], 'Initializing parallelization FVA subproblem tasks on '+str(nOfProcessors)+' processors\n')
        try:
            # Parallelize
            size = len(fluxes)
            fvaData = []
            p = Pool(nOfProcessors)
            func = partial(fvaSubProblem, Oinstance, OrOpt, size, optV, lowerBound, upperBound)
            r = p.map_async(func, fluxes, callback=fvaData.append)
            r.wait()
        except:
            mesg(settings['logFile'], 'Error: could not parallelize FVA subproblems\n')
            return fvaData

        # Add the model and its settings to fvaData
        try:
            fvaData[0].append(model + ' ' + str(settings))
        except:
            mesg(settings['logFile'], 'Error: could not append settings to FVA data list. Check if output is empty')

        if fileName != '':
            try:
                # saves fva using pickle
                with open(fileName, 'wb') as f:
                    pickle.dump(fvaData, f)
            except:
                mesg(settings['logFile'], 'Error: could not save FVA data list to file')
                return fvaData

        # print where FVA is saved and time duration
        mesg(settings['logFile'], '\nDone. FVA results saved as ' + fileName + ' using the pickle package. Duration of FVA: ' + str(round(time.time() - start_time))+' s')
        return fvaData
