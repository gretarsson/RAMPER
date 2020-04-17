from __future__ import division  # safety with double division
from scipy.io import loadmat, savemat
from scipy import sparse, arange, unique
import numpy
import os
from pyomo.environ import *
from pyomo.opt import SolverFactory, SolverStatus, TerminationCondition
import logging
import cobra
import time
import pickle  # debugging
# TODO: 1. make assumptions clear: prot_ nomenclature; all prot_ pseudo-metabolites and pseudoreactions must be last columnwise and rowwise; prot_pool as last row AND last column
#       2. Find a way to find the b vector quickly when the input is of cobrapy type. As of now it is assumed that the b vector is a zero vector

#
# Returns a list of None
#
def declareNone(numNullReturns):
    noneList = [None for i in range(numNullReturns)]
    return tuple(noneList)


def mesg(logFile, mesgToAdd):
    fl = open(logFile, "a+")
    fl.write(mesgToAdd)
    fl.close()

# noinspection PyBroadException
def parseData(matFile, settings):
    mesg(settings['logFile'], 'In parseData\n')
    # Initiation
    Sspr, Kspr, P, MW, bS, bK, cS, ubS, lbS, ubK, lbK, totMass, slmInd, Kmap, Smap, SmodMap = declareNone(16)
    stat = 1  # assume failure
    isMat = False  # flag for input
    isCobra = False  # flag for input

    # determine input type -> .xml\cobrapy or .mat
    try:
        if str(matFile)[-4:] == '.mat':
            isMat = True
            mesg(settings['logFile'], '\tModel input determined to a matlab file.\n')
            pathToFile = settings['pathToFile']
            matModName = settings['matModName']
            matFile = pathToFile + matFile
            # load GECKO data from mat file
            if os.path.exists(matFile):
                modInfo = loadmat(matFile)
                ## delete temporary .mat file created in the case where a cobrapy model object is the input
                #if convFlag:
                    #os.remove(matFile)
                mesg(settings['logFile'], '\tDone reading ' + matFile + '\n')
                # Look for matModName if not given as argument
                if not matModName:  # TODO: what if model does not have __version__ and so forth?
                    matModName = \
                        [modInfo.keys()[i] for i in range(len(modInfo.keys())) if modInfo.keys()[i] != '__version__' and
                         modInfo.keys()[i] != '__header__' and modInfo.keys()[i] != '__globals__'][0]
                    if not matModName:
                        mesg(settings['logFile'], '\tError: Could not find name of model in matlab structure\n')
                        return Sspr, Kspr, P, MW, bS, bK, cS, ubS, lbS, ubK, lbK, totMass, slmInd, Kmap, Smap, SmodMap, stat
                else:
                    if not (matModName in modInfo):
                        mesg(settings['logFile'],
                             '\tError: Could not find \'' + str(matModName) + '\' as key in model structure\n')
                        return Sspr, Kspr, P, MW, bS, bK, cS, ubS, lbS, ubK, lbK, totMass, slmInd, Kmap, Smap, SmodMap, stat
            else:
                mesg(settings['logFile'], '\tError: Matlab File ' + matFile + ' does not exist\n')
                return Sspr, Kspr, P, MW, bS, bK, cS, ubS, lbS, ubK, lbK, totMass, slmInd, Kmap, Smap, SmodMap, stat
            if set(['S', 'b', 'c', 'ub', 'lb', 'metNames', 'rxnNames']).issubset(set(modInfo[matModName].dtype.names)):
                Smod = modInfo[matModName]['S']  # modified S mat with kcat submat
                Smod = sparse.csr_matrix(Smod[0][0])
                bmod = modInfo[matModName]['b']  # modified b
                cmod = modInfo[matModName]['c']  # modified c
                ubmod = modInfo[matModName]['ub']  # modified ub
                lbmod = modInfo[matModName]['lb']  # modeified lb
                metName = modInfo[matModName]['metNames']  # list of metabolite ( + protein) names
                rxnName = modInfo[matModName]['rxnNames']  # list of reaction names
                mesg(settings['logFile'], '\tDone extracting data LP data from ' + matFile + '\n')
            else:
                mesg(settings['logFile'], '\tError: Mat structure ' + matModName + ' fails to have all structural elements\n')
                return Sspr, Kspr, P, MW, bS, bK, cS, ubS, lbS, ubK, lbK, totMass, Kmap, Smap, SmodMap, stat

                #
            # Put this in a more standard form for us
            #
            # The i,j - element of the modified S matrix is Smod[0][0][i][j]
            S = []  # non-kcat submatrix
            K = []  # kcat submatrix
            bS = []  # non-kcat subvector
            bK = []  # right-hand side for K matrix
            metS = []  # met names for the stoichiometric system
            metK = []  # met names for the enzyme system
            poolRow = 0  # 0 - no pool constr, 1 - pool constr
            
	    # initiates uncertainty matrix
            uncMat = settings['uncertaintyMatrix']
	    uncMat = numpy.array(uncMat)
	    uncMatBool = uncMat.any()
            # maps the row-indeces from K and S to the uncertainty matrix
 	    Kmap = []
	    Smap = []
	    SmodMap = []
	    Sind = 0
	    Kind = 0
            try:
                # Segment the rows into stoichiometry and enzyme
                for i in range(Smod.shape[0]):
                    if metName[0][0][i][0][0].find('prot_') < 0:
                        # S.append(Smod[0][0][i])
                        S.append(Smod[i, :].toarray()[0])
                        bS.append(bmod[0][0][i][0])
                        metS.append(metName[0][0][i][0][0])
			# append index mappings
			SmodMap.append(['S',Sind])
			Sind += 1
			if uncMatBool:
				Smap.append(i)
                    else:
                        # K.append(Smod[0][0][i])
                        K.append(Smod[i, :].toarray()[0])
                        bK.append(bmod[0][0][i][0])
                        metK.append(metName[0][0][i][0][0])
                        if metName[0][0][i][0][0].find('prot_pool') >= 0:  # this is pools line
                            poolRow = 1
                            mesg(settings['logFile'], '\tPool constraint found\n')
			# append index mappings
			SmodMap.append(['K',Kind])
			Kind += 1
			if uncMatBool:
				Kmap.append(i)
                # find the vertical column sep between flux and enzyme
                # colSepInd = [j for j in range(len(Smod[0][0][poolRowInd])) if Smod[0][0][poolRowInd][j] != 0][0]
                colSepInd = Smod.shape[1] - len(K)

                # Compress S and K to sparse matrices
                Sspr = sparse.csr_matrix(S)
                Kspr = sparse.csr_matrix(K)
                mesg(settings['logFile'], '\tDone parsing the larger S and K from the GECKO LP matrix\n')
		if uncMatBool:
			mesg(settings['logFile'], '\tDone creating index mapping between the uncertainty matrix and the parsed S and K matrices\n')
            except:
                mesg(settings['logFile'], '\tError: unable to identify S and K from mat file\n')
                return Sspr, Kspr, P, MW, bS, bK, cS, ubS, lbS, ubK, lbK, totMass, slmInd, Kmap, Smap, SmodMap, stat
        # check if input is cobrapy model
        elif type(matFile) is cobra.core.model.Model:
            isCobra = True
            mesg(settings['logFile'], '\tModel input determined to be cobrapy model object\n')
            try:
                # extract S matrix
                Smod = cobra.util.create_stoichiometric_matrix(matFile)
                Smod = sparse.csr_matrix(Smod)
                mesg(settings['logFile'], '\tDone extracting S matrix\n')

                # adjust nonexisting--listed as NoneType--variable bounds to numpy infinities as NoneType is not handled by cobra.util.constraint_matrices.
                # also note that setting non-inf bounds may cause solver instability.
                for i in range(len(matFile.variables)//2):
                    if matFile.variables[2*i].ub is None:
                        matFile.variables[2*i].ub = numpy.Inf
                    if matFile.variables[2*i].lb is None:
                        matFile.variables[2*i].lb = numpy.NINF
                mesg(settings['logFile'], '\tDone updating cobrapy variable bounds\n')

                # extract b vector
                #conMats = cobra.util.constraint_matrices(matFile)  # this takes a long time... we only really use the b
                mesg(settings['logFile'], '\tConstraint matrices generated\n')
                #bmod = [conMats.b[j] for j in range(len(conMats.b))]
                bmod = [0]*len(matFile.metabolites)
                mesg(settings['logFile'], '\tb vector found\n')

                # extract the bound vectors
                # the variables are split up into irreversible reactions by cobrapy (even if they are already split)
                ubmod = []
                lbmod = []
                for i in range(len(matFile.variables)//2):
                    ubmod.append(matFile.variables[2*i].ub)  # upper bound of var i
                    if matFile.variables[2*i].lb == 0:  # if lower bound is positive we have to ignore the lower bound of the reversed reaction
                        lbmod.append(-matFile.variables[2 * i + 1].ub)
                    else:
                        lbmod.append(matFile.variables[2*i].lb)
                mesg(settings['logFile'], '\tVar bound vectors found\n')
                #ubmod = [conMats.variable_bounds[2*j][1] for j in range(len(conMats.variable_bounds)//2)]
                #lbmod = [conMats.variable_bounds[2*j][0] for j in range(len(conMats.variable_bounds)//2)]
                mesg(settings['logFile'], '\tDone extracting b vector and variable bound vectors\n')

                # generate c vector
                cmod = []
                objDict = matFile.objective.get_linear_coefficients(matFile.variables)
                for i in range(len(matFile.variables)//2):
                    cmod.append(objDict[matFile.variables[2*i]])
                mesg(settings['logFile'], '\tDone extracting c vector\n')

                # get metabolite names
                metName = [matFile.metabolites[i].name for i in range(len(matFile.metabolites))]
                mesg(settings['logFile'], '\tDone extracting metabolite names\n')

                # get reaction names
                rxnName = [matFile.reactions[i].name for i in range(len(matFile.reactions))]
                mesg(settings['logFile'], '\tDone extracting reaction names\n')
                mesg(settings['logFile'], '\tDone extracting data from ' + str(matFile) + ' cobrapy model\n')
            except:
                mesg(settings['logFile'], '\tError: Could not extract data successfully\n')
                return Sspr, Kspr, P, MW, bS, bK, cS, ubS, lbS, ubK, lbK, totMass, slmInd, Kmap, Smap, SmodMap,stat

            # Put this in a more standard form for us
            #
            # The i,j - element of the modified S matrix is Smod[0][0][i][j]
            S = []  # non-kcat submatrix
            K = []  # kcat submatrix
            bS = []  # non-kcat subvector
            bK = []  # right-hand side for K matrix
            metS = []  # met names for the stoichiometric system
            metK = []  # met names for the enzyme system
            poolRow = 0  # 0 - no pool constr, 1 - pool constr

	    # initiates uncertainty matrix
            uncMat = settings['uncertaintyMatrix']
	    uncMat = numpy.array(uncMat)
	    uncMatBool = numpy.any(uncMat)
            # maps the row-indeces from K and S to the uncertainty matrix
 	    Kmap = []
	    Smap = []
	    SmodMap = []
	    Sind = 0
	    Kind = 0
		
            try:
                # Segment the rows into stoichiometry and enzyme
                for i in range(Smod.shape[0]):
                    if metName[i].find('prot_') < 0:
                        # S.append(Smod[0][0][i])
                        S.append(Smod[i, :].toarray()[0])
                        bS.append(bmod[i])
                        metS.append(metName[i])
			# append index mappings
			SmodMap.append(['S',Sind])
			Sind += 1
			if uncMatBool:
				Smap.append(i)
                    else:
                        # K.append(Smod[0][0][i])
                        K.append(Smod[i, :].toarray()[0])
                        bK.append(bmod[i])
                        metK.append(metName[i])
                        if metName[i].find('prot_pool') >= 0:  # this is pools line
                            poolRow = 1
                            mesg(settings['logFile'], '\tPool constraint found\n')
			# append index mappings
			SmodMap.append(['K',Kind])
			Kind += 1
			if uncMatBool:
				Kmap.append(i)
                # find the vertical column sep between flux and enzyme
                # colSepInd = [j for j in range(len(Smod[0][0][poolRowInd])) if Smod[0][0][poolRowInd][j] != 0][0]
                colSepInd = Smod.shape[1] - len(K)

                # Compress S and K to sparse matrices
                Sspr = sparse.csr_matrix(S)
                Kspr = sparse.csr_matrix(K)
                mesg(settings['logFile'], '\tDone parsing the larger S and K from the GECKO LP matrix\n')
		if uncMatBool:
			mesg(settings['logFile'], '\tDone creating index mapping between the uncertainty matrix and the parsed S and K matrices\n')
            except:
                mesg(settings['logFile'], '\tError: unable to identify S and K from mat file\n')

        else:
            mesg(settings['logFile'], '\tError: Could not determine model input type')
            return Sspr, Kspr, P, MW, bS, bK, cS, ubS, lbS, ubK, lbK, totMass, slmInd, Kmap, Smap, SmodMap, stat
    except:
        mesg(settings['logFile'], '\tError: Could not determine input type of model\n')
        return Sspr, Kspr, P, MW, bS, bK, cS, ubS, lbS, ubK, lbK, totMass, slmInd, Kmap, Smap, SmodMap, stat

    try:
        # Snatch the submatrices for our model formulation
        Sspr = Sspr[0:Sspr.shape[0], 0:colSepInd]
        P = Kspr[0:Kspr.shape[0] - poolRow, colSepInd:Kspr.shape[1] - poolRow]
        P = P * arange(1., P.shape[1] + 1)  # return a permutation vector
        if len(unique(P)) != len(P):
            mesg(settings['logFile'], '\tError: Permutation test failed\n')
            return Sspr, Kspr, P, MW, bS, bK, cS, ubS, lbS, ubK, lbK, totMass, slmInd, Kmap, Smap, SmodMap, stat
        if poolRow:
            MW = numpy.absolute([Kspr[Kspr.shape[0] - 1, colSepInd:Kspr.shape[1] - 1].toarray()][0][0])  # absolute value in case coefficients are neg.
        else:
            MW = [0] * Kspr.shape[0]
        Kspr = numpy.absolute(Kspr[0:Kspr.shape[0] - poolRow, 0:colSepInd])  # absolute value to get positive kcat vals
        mesg(settings['logFile'], '\tDone parsing S, K, P, and MW\n')
    except:
        mesg(settings['logFile'], '\tError: unable to parse S, K, MW, and P from sparse matrix\n')
        return Sspr, Kspr, P, MW, bS, bK, cS, ubS, lbS, ubK, lbK, totMass, slmInd, Kmap, Smap, SmodMap, stat

    if isMat:
        try:
            # Finish off by identifying and grabbing the column data
            Srange = [j for j in range(colSepInd)]
            Krange = [j for j in range(len(cmod[0][0]) - poolRow) if j >= colSepInd]
            cS = [cmod[0][0][j][0] for j in Srange]
            ubS = [ubmod[0][0][j][0] for j in Srange]
            lbS = [lbmod[0][0][j][0] for j in Srange]
            ubK = [ubmod[0][0][j][0] for j in Krange]
            lbK = [lbmod[0][0][j][0] for j in Krange]
            if poolRow:
                totMass = ubmod[0][0][len(cmod[0][0]) - 1][0]
            else:
                totMass = 0
            mesg(settings['logFile'], '\tDone parsing cS, ubS, lbS, ubK, lbK, and totMass\n')
        except:
            mesg(settings['logFile'], '\tError: unable to parse column information\n')
            return Sspr, Kspr, P, MW, bS, bK, cS, ubS, lbS, ubK, lbK, totMass, slmInd, Kmap, Smap, SmodMap, stat

        try:
            # Find SLIME reaction/column indices
            #slmInd = [j for j in range(colSepInd) if 'SLIME' in rxnName[0][0][j][0][0]]
            slmInd = []
            for j in range(colSepInd):
                # avoids error message if reaction name is empty and prints to logFile that this is the case
                if rxnName[0][0][j][0].size != 0:
                    strRxnName = rxnName[0][0][j][0][0]
                else:
                    strRxnName = ''
                    mesg(settings['logFile'], '\tReaction at index ' + str(j) + ' has no name\n')
                if 'SLIME' in strRxnName:
                    slmInd.append(j)
            mesg(settings['logFile'], '\tDone identifying SLIME reaction columns. ' + str(len(slmInd)) + ' SLIME reactions found\n')
        except:
            mesg(settings['logFile'], '\tError: unable to identify SLIME reaction columns\n')
            return Sspr, Kspr, P, MW, bS, bK, cS, ubS, lbS, ubK, lbK, totMass, slmInd, Kmap, Smap, SmodMap, stat
    elif isCobra:
        try:
            # Finish off by identifying and grabbing the column data
            Srange = [j for j in range(colSepInd)]
            Krange = [j for j in range(len(cmod) - poolRow) if j >= colSepInd]
            cS = [cmod[j] for j in Srange]
            ubS = [ubmod[j] for j in Srange]
            lbS = [lbmod[j] for j in Srange]
            ubK = [ubmod[j] for j in Krange]
            lbK = [lbmod[j] for j in Krange]
            if poolRow:
                totMass = ubmod[len(cmod) - 1]
            else:
                totMass = 0
            mesg(settings['logFile'], '\tDone parsing cS, ubS, lbS, ubK, lbK, and totMass\n')
        except:
            mesg(settings['logFile'], '\tError: unable to parse column information\n')
            return Sspr, Kspr, P, MW, bS, bK, cS, ubS, lbS, ubK, lbK, totMass, slmInd, Kmap, Smap, SmodMap, stat

        try:
            # Find SLIME reaction/column indices
            # slmInd = [j for j in range(colSepInd) if 'SLIME' in rxnName[0][0][j][0][0]]
            slmInd = []
            for j in range(colSepInd):
                # avoids error message if reaction name is empty and prints to logFile that this is the case
                if len(rxnName[j]) != 0:
                    strRxnName = rxnName[j]
                else:
                    strRxnName = ''
                    mesg(settings['logFile'], '\tReaction at index ' + str(j) + ' has no name\n')
                if 'SLIME' in strRxnName:
                    slmInd.append(j)
            mesg(settings['logFile'], '\tDone identifying SLIME reaction columns. ' + str(len(slmInd)) + ' SLIME reactions found\n')
        except:
            mesg(settings['logFile'], '\tError: unable to identify SLIME reaction columns\n')
            return Sspr, Kspr, P, MW, bS, bK, cS, ubS, lbS, ubK, lbK, totMass, slmInd, Kmap, Smap, SmodMap, stat

    # save variables (debugging)
    #with open('pyObj.pkl', 'w') as f:
        #pickle.dump({'Sspr': Sspr, 'Kspr': Kspr, 'bS': bS, 'bK': bK, 'ubS': ubS, 'lbS': lbS, 'ubK': ubK, 'lbK': lbK, 'cS': cS, 'colSepInd': colSepInd, 'P': P, 'MW': MW, 'totMass': totMass, 'slmInd': slmInd}, f)

    # clear some memory
    del Smod, bmod, cmod, ubmod, lbmod, S, K, metName
    if isMat:
        del modInfo

    mesg(settings['logFile'], '\tParsing complete\n')
    stat = 0  # we have not failed
    return Sspr, Kspr, P, MW, bS, bK, cS, ubS, lbS, ubK, lbK, totMass, slmInd, Kmap, Smap, SmodMap, stat


def modelKcatUncertainty(Kspr, Kmap, settings):
    mesg(settings['logFile'], 'In modelKcatUncertainty\n')

    # Initialize
    Rspr, Rindx = declareNone(2)
    stat = 1
    numRows = 0
    Rindx = [0]
    uniqueVal = {}
    nOfUncPar = [0 for i in range(Kspr.shape[0])]
    linK = settings['linearKUncertainty']
    uncMat = numpy.array(settings['uncertaintyMatrix'])
    uncMatBool = uncMat.any()
    covars = settings['covarianceList']

    try:
        # first check if model had a K matrix in the first place
        if min(Kspr.shape[0], Kspr.shape[1]) != 0:
            try:
                if linK:  # TODO: this whole thing if only uncmat
                    minKUnc = settings['minKUnc']
                    maxKUnc = settings['maxKUnc']
                    percKForced = settings['percKForced']
                    # find upper and lower bounds on kcats
                    minK = 1 / Kspr.max()
                    maxK = 1 / Kspr.data.min()
		if uncMatBool:
			pass
                else:
                    # Verify (and lengthen) maxPercUncK list TODO: make conditional statement nicer
                    maxPercUncK = settings['maxPercUncK']
                    if type(maxPercUncK) is not list or type(maxPercUncK) is list and len(maxPercUncK) != 1 \
                            and len(maxPercUncK) != Kspr.shape[0]:
                        mesg(settings['logFile'], '\tError: argument maxPercUncK must be list of size 1 or ' + str(Kspr.shape[0])
                             + ' (n of rows in K) \n')
                        return Rspr, Rindx, stat
                    if len(maxPercUncK) == 1:
                        maxPercUncK = [maxPercUncK[0]] * Kspr.shape[0]
            except:
                mesg(settings['logFile'], '\tError: Failed checking arguments\n')
                return Rspr, Rindx, stat

            try:
                # first pass of rows of K to index R and to calculate its size
                for i in range(Kspr.shape[0]):
		    if uncMatBool:
			Mi = Kmap[i]  # switch to index in Smod
		    	nOfUncPar[i] = numpy.count_nonzero(uncMat[Mi])
			if covars:
				covarRow = covars[Mi]
				if covarRow:
					for l in range(len(covarRow)):
						nOfUncPar[i] = nOfUncPar[i] - (len(covarRow[l]) - 1)
		    else:
                   	uniqueVal[i] = unique(Kspr[i, :].toarray())
                   	uniqueVal[i] = [uniqueVal[i][j] for j in range(len(uniqueVal[i])) if uniqueVal[i][j] != 0]
			nOfUncPar[i] = len(uniqueVal[i])
		    numRows = numRows + nOfUncPar[i]
		    Rindx.append(Rindx[len(Rindx) - 1] + nOfUncPar[i])
                mesg(settings['logFile'], '\tCompleted first pass through K\n')
            except:
                mesg(settings['logFile'], '\tError: failed in first pass through K\n')
                return Rspr, Rindx, stat

            try:
                # initialize R
                R = numpy.zeros((numRows, Kspr.shape[1]))
                k = 0
                # build R
                for i in range(Kspr.shape[0]):
		    # pipeline for uncertainty matrix input	
		    if uncMatBool:
			uncRow = uncMat[Kmap[i], 0:Kspr.shape[1]] # as Kspr is cut short columnwise
			if covars:
				covar = covars[Kmap[i]]
			else:
				covar = []
			uncInd = numpy.nonzero(uncRow)[0]
			tempR = numpy.array( [[0 for col in range(Kspr.shape[1])] for row in range(len(uncInd))] ).astype(float)
			tempk = 0 
			# skip if no uncertain parameters in row
			if uncInd.size:
				# first treat all uncorrelated then we merge the R afterwards
				for j in uncInd: 
					if uncRow[j] == -1:
						
						mesg(settings['logFile'], 'd')
						if Kspr[i,:].toarray()[i][j] == 0:
							mesg(settings['logFile'], '\tZero coefficient in K assigned -1. Division by zero ensues.\n')	
						kcat = 1/Kspr[i,:].toarray()[i][j]
						percK = (kcat - minK) / (maxK-minK)*percKForced[i]
						uncK = minKUnc + (maxKUnc - minKUnc) * percK
						tempR[tempk, j] = uncK / ((1 + uncK) * kcat)
					else:
						mesg(settings['logFile'], 'tempR = '+str(tempR[tempk, j])+'\tuncRow[j] = '+str(uncRow[j])+'\n')
						tempR[tempk, j] = uncRow[j]
					tempk = tempk + 1
								
				# find correlated rows
				rowsMergeL = [[] for row in range(len(covar))]
				for j in range(len(covar)):
					for s in range(len(covar[j])):
						rowInd = numpy.nonzero(tempR[:,covar[j][s]])[0][0]
						rowsMergeL[j].append(rowInd)	
				# merge them
				for r in range(len(rowsMergeL)):
					rowsMerge = numpy.array(rowsMergeL[r])
					masterRow = numpy.sum(tempR[rowsMerge,:],axis=0)		
					tempR = numpy.delete(tempR, rowsMerge,axis=0)
					tempR = numpy.vstack((tempR, masterRow))
				# add the tempR to R
				tempRlen = tempR.shape[0]
				R[k:(k+tempRlen),:] = tempR		
				k = k + tempRlen
		    else:
			# pipeline without uncertainty matrix input
                    	if not linK:
                        	maxPercUncKi = maxPercUncK[i]
                    	for j in range(len(uniqueVal[i])):	
                        	tmpInd = numpy.nonzero(Kspr[i, :] == uniqueVal[i][j])[1]
                        	if not linK:
                            		R[k, tmpInd] = maxPercUncKi * uniqueVal[i][j]
                        	else:
                            		kcat = 1/uniqueVal[i][j]
                            		percK = (kcat - minK) / (maxK-minK)*percKForced[i]
                            		uncK = minKUnc + (maxKUnc - minKUnc) * percK
                            		R[k, tmpInd] = uncK / ((1 + uncK) * kcat)
                        	k = k + 1
                mesg(settings['logFile'], '\tCompleted construction of RK\n')
		mesg(settings['logFile'], '\tR final :\n\t'+str(R)+'\n')
            except:
                mesg(settings['logFile'], '\tError: failed in second pass through K\n')
                return Rspr, Rindx, stat

            # Reduce storage
            Rspr = sparse.csr_matrix(R)
            del R, uniqueVal

            # we are done
            mesg(settings['logFile'], '\tUncertainty model complete\n')
            stat = 0  # we did not fail
            return Rspr, Rindx, stat
        else:
            # No K matrix found
            mesg(settings['logFile'], '\tDid not find K matrix although K uncertainty was set to be modelled. K uncertainty modelling was skipped\n')
            stat = 0  # we did not fail
            return Rspr, Rindx, stat
    except:
        mesg(settings['logFile'], 'K matrix has no shape() method.')
        stat = 1  # we did fail
        return Rspr, Rindx, stat


def modelSUncertainty(Sspr, bS, slmInd, Smap,  settings):
    mesg(settings['logFile'], 'In modelSUncertainty\n')

    # Initialize returns
    SsprCer, SsprUnc, bSCer, bSUnc, Rspr, Rindx = declareNone(6)
    stat = 1  # assume failure
    maxPercUncS = settings['maxPercUncS']
    uncMat = numpy.array(settings['uncertaintyMatrix'])
    uncMatBool = uncMat.any()
    covars = settings['covarianceList']

    # initialize some structures
    cerRowInd = []
    uncerRowInd = []
    Rindx = [0]
    nOfUncPar = {}
    numRows = 0
    try:
        # first pass of rows of S to calculate (un)certain indices 
        for i in range(Sspr.shape[0]):
	    if uncMatBool: 
			Mi = Smap[i]  # switch to index in Smod
		    	nOfUncPar[i] = numpy.nonzero(uncMat[Mi])[0]
	    else:
		    if not slmInd:  # checks if there is no SLIME reactions
			nOfUncPar[i] = unique(Sspr[i, :].toarray()[0])
		    else:
			nOfUncPar[i] = Sspr[i, :].toarray()[0]
			nOfUncPar[i][numpy.array(slmInd)] = 0  # treat SLIME columns as 0, i.e. certain, slmInd cannot be []
			nOfUncPar[i] = unique(nOfUncPar[i])
		    nOfUncPar[i] = [nOfUncPar[i][j] for j in range(len(nOfUncPar[i])) if (nOfUncPar[i][j] % 1 != 0) and (nOfUncPar[i][j] % 1 != 0.5)]
            if len(nOfUncPar[i]) == 0:  # This is a certain row
                cerRowInd.append(i)
            else:
                uncerRowInd.append(i)
                numRows = numRows + len(nOfUncPar[i])
		if covars:
			covarRow = covars[Mi]
			if covarRow:
				mesg(settings['logFile'], '\tCovariance found in S matrix. May lead to error in relaxed-SS constraint implementation\n')
				for l in range(len(covarRow)):
					numRows = numRows - (len(covarRow[l]) - 1)
                Rindx.append(Rindx[len(Rindx) - 1] + len(nOfUncPar[i]))
        mesg(settings['logFile'], '\tCompleted first pass through S. ' + str(numRows) + ' uncertain rows created\n')
    except:
        mesg(settings['logFile'], '\tError: Failed in first pass through S\n')
        return SsprCer, SsprUnc, bSCer, bSUnc, Rspr, Rindx, stat

    try:
        # Divide S and b into their certian and uncertain parts
        SsprCer = Sspr[cerRowInd, :]
        SsprUnc = Sspr[uncerRowInd, :]
        bSCer = [bS[cerRowInd[j]] for j in range(len(cerRowInd))]
        bSUnc = [bS[uncerRowInd[j]] for j in range(len(uncerRowInd))]
        mesg(settings['logFile'], '\tParsed S into (un)certain parts\n')
    except:
        mesg(settings['logFile'], '\tError: Failed to parse S into (un)certain parts\n')
        return SsprCer, SsprUnc, bSCer, bSUnc, Rspr, Rindx, stat

    try:
        # build R
        R = numpy.zeros((numRows, Sspr.shape[1]))
        k = 0
        for i in range(len(uncerRowInd)):
	    # pipeline for uncertainty matrix input	
	    if uncMatBool:
		uncRow = uncMat[Smap[uncerRowInd[i]]]
		if covars:
			covar = covars[Smap[i]]
		else:
			covar = []
		uncInd = numpy.nonzero(uncRow)[0]
		tempR = numpy.array( [[0 for col in range(Sspr.shape[1])] for row in range(len(uncInd))] ).astype(float)
		tempk = 0 
		# skip if no uncertain parameters in row
		if uncInd.size:
			# first treat all uncorrelated then we merge the R afterwards
			for j in uncInd: 
				if uncRow[j] == -1:
					tempR[tempk, j] = maxPercUncS
				else:
					tempR[tempk, j] = uncRow[j]
				tempk = tempk + 1
			# find correlated rows
			rowsMergeL = [[] for row in range(len(covar))]
			for j in range(len(covar)):
				for s in range(len(covar[j])):
					rowInd = numpy.nonzero(tempR[:,covar[j][s]])[0][0]
					rowsMergeL[j].append(rowInd)	
			# merge them
			for r in range(len(rowsMergeL)):
				rowsMerge = numpy.array(rowsMergeL[r])
				masterRow = numpy.sum(tempR[rowsMerge,:],axis=0)		
				tempR = numpy.delete(tempR, rowsMerge,axis=0)
				tempR = numpy.vstack((tempR, masterRow))
			# add the tempR to R
			tempRlen = tempR.shape[0]
			R[k:(k+tempRlen),:] = tempR		
			k = k + tempRlen
	    else:
		    for j in range(len(nOfUncPar[uncerRowInd[i]])):
			tmpInd = numpy.nonzero(Sspr[uncerRowInd[i], :] == nOfUncPar[uncerRowInd[i]][j])[1]
			R[k, tmpInd] = maxPercUncS * nOfUncPar[uncerRowInd[i]][j]
			k = k + 1
        mesg(settings['logFile'], '\tCompleted construction of RS\n')
    except:
        mesg(settings['logFile'], '\tError: Failed to build RS\n')
        return SsprCer, SsprUnc, bSCer, bSUnc, Rspr, Rindx, stat
    # Reduce storage
    Rspr = sparse.csr_matrix(R)
    del R, nOfUncPar, #cerRowInd, uncerRowInd
    # we are done
    mesg(settings['logFile'], '\tUncertainty Model Complete\n')
    stat = 0  # we succeded
    return SsprCer, SsprUnc, bSCer, bSUnc, Rspr, Rindx, cerRowInd, uncerRowInd, stat


def solveModel(SsprCer, SsprUnc, Kspr, P, MW, bSCer, bSUnc, \
               bK, cS, ubS, lbS, ubK, lbK, totMass, \
               RSspr, RSindx, RKspr, RKindx, M, cerRowInd, uncerRowInd, settings):
    mesg(settings['logFile'], 'In solveModel\n')

    # Initiate
    sol = declareNone(1)
    stat = 1  # assume failure

    if not settings['seeWarnings']:  # Turn off Warnings - pesky little things
        logging.getLogger('pyomo.core').setLevel(logging.ERROR)

    try:
        # Set some flags for types of uncertainty we are dealing with
        #   these assume certainty is distinguished by R being [] as the function is called
        if sparse.issparse(RSspr):
            SUncFlag = (RSspr.shape[0] != 0)
        else:
            if RSspr is None:
                SUncFlag = False
            else:
                SUncFlag = (len(RSspr) != 0)
        if sparse.issparse(RKspr):
            KUncFlag = (RKspr.shape[0] != 0)
        else:
            if RKspr is None:
                KUncFlag = False
            else:
                KUncFlag = (len(RKspr) != 0)
        mesg(settings['logFile'], '\tIdentified uncertainty structure via RS and Rk\n')
    except:
        mesg(settings['logFile'], '\tError: unable to discern uncertainty structure via RS and Rk\n')
        return sol, stat

    try:
        # idenfity nonzero elements of S and K to expidite constraint generation
        nonZeroIndSCer = numpy.nonzero(SsprCer != 0)  # Identify the nonzero elements SsprCer
        nonZeroIndK = numpy.nonzero(Kspr != 0)  # Identify the nonzero elements
        if SUncFlag:
            nonZeroIndSUnc = numpy.nonzero(SsprUnc != 0)  # Identify the nonzero elements SsprUnc
            nonZeroIndRS = numpy.nonzero(RSspr != 0)  # Identify the nonzero elements RSspr
        if KUncFlag:
            nonZeroIndRK = numpy.nonzero(RKspr != 0)  # Identify the nonzero elements
        mesg(settings['logFile'], '\tIdentified nonzero elements of S and K\n')
    except:
        mesg(settings['logFile'], '\tError: Failed to identified nonzero elements of S and K\n')
        return sol, stat

    #
    # Define the model
    #
    # model def
    yst = AbstractModel()

    try:
        # Sets
        yst.metSCer = RangeSet(0, SsprCer.shape[0] - 1)  # index set for rows of SCer
        yst.flux = RangeSet(0, SsprCer.shape[1] - 1)  # index set for cols of SCer (number of fluxes)
        # checks if there is any enzymes
        if min(Kspr.shape[0], Kspr.shape[1]) > 0:
            yst.enzK = RangeSet(0, Kspr.shape[0] - 1)  # index set for rows of K
        if SUncFlag:
            yst.metSUnc = RangeSet(0, SsprUnc.shape[0] - 1)  # index set for rows of SUnc
            yst.RSxid = RangeSet(0, RSspr.shape[0] - 1)  # index set for RSx vars
        if KUncFlag:
            yst.RKxid = RangeSet(0, RKspr.shape[0] - 1)  # index set for Rx vars
        mesg(settings['logFile'], '\tDone instantiating sets\n')
    except:
        mesg(settings['logFile'], '\tError: Failed to instantiate sets\n')
        return sol, stat

    try:
        # Variables
        def ensureBoundsS(yst, j):
            return (lbS[j], ubS[j])

        yst.v = Var(yst.flux, within=Reals, bounds=ensureBoundsS)  # flux variables

        def ensureBoundsK(yst, j):
            return (lbK[j], ubK[j])

        # Check if there is any enzymes
        if min(Kspr.shape[0], Kspr.shape[1]) > 0:
            yst.e = Var(yst.enzK, within=Reals, bounds=ensureBoundsK)  # enzyme variables
        if SUncFlag:
            yst.Su = Var(yst.RSxid, within=Reals)  # Rx variables
            yst.Sw = Var(yst.metSUnc, within=NonNegativeReals)  # bound variables for ||RSUnc x||
        if KUncFlag:
            yst.Ku = Var(yst.RKxid, within=Reals)  # Rx variables
            yst.Kw = Var(yst.enzK, within=NonNegativeReals)  # bound variables for ||RKUnc x||
        mesg(settings['logFile'], '\tDone instantiating variables\n')
    except:
        mesg(settings['logFile'], '\tError: Failed to instantiate variables\n')
        return sol, stat

    try:
        # Objective
        def calcGrowth(yst):
            return sum(cS[j] * yst.v[j] for j in yst.flux)

        yst.grwthRate = Objective(rule=calcGrowth, sense=maximize)
        mesg(settings['logFile'], '\tDefined the objective \n')
    except:
        mesg(settings['logFile'], '\tError: Failed to define the objective\n')
        return sol, stat

    try:
        # Steady state constraints for certain rows
        def ensureSteadyStateCer(yst, i):
            rowInd = numpy.nonzero(nonZeroIndSCer[0] == i)[0]
            return sum(SsprCer[i, nonZeroIndSCer[1][k]] * yst.v[nonZeroIndSCer[1][k]] for k in rowInd) == bSCer[i]

        yst.steadyStateCer = Constraint(yst.metSCer, rule=ensureSteadyStateCer)

        # Uncertain steady state constraints
        if SUncFlag:
            def ensureSteadyStateUncUpper(yst, i):
                rowInd = numpy.nonzero(nonZeroIndSUnc[0] == i)[0]
                return sum(SsprUnc[i, nonZeroIndSUnc[1][k]] * yst.v[nonZeroIndSUnc[1][k]] for k in rowInd) \
                       + yst.Sw[i] <= bSUnc[i] + M[i]

            yst.steadyStateUncUpper = Constraint(yst.metSUnc, rule=ensureSteadyStateUncUpper)

            def ensureSteadyStateUncLower(yst, i):
                rowInd = numpy.nonzero(nonZeroIndSUnc[0] == i)[0]
                return sum(SsprUnc[i, nonZeroIndSUnc[1][k]] * yst.v[nonZeroIndSUnc[1][k]] for k in rowInd) \
                       - yst.Sw[i] >= -bSUnc[i] - M[i]

            yst.steadyStateUncLower = Constraint(yst.metSUnc, rule=ensureSteadyStateUncLower)

            def equateRobustSUnc(yst, i):
                rowInd = numpy.nonzero(nonZeroIndRS[0] == i)[0]
                return sum(RSspr[i, nonZeroIndRS[1][k]] * yst.v[nonZeroIndRS[1][k]] for k in rowInd) - yst.Su[i] == 0

            yst.rbsBalSUnc = Constraint(yst.RSxid, rule=equateRobustSUnc)

            def enforceRobustSUnc(yst, i):
                return sum(yst.Su[k] * yst.Su[k] for k in range(RSindx[i], RSindx[i + 1])) - yst.Sw[i] * yst.Sw[i] <= 0

            yst.socpConstSUnc = Constraint(yst.metSUnc, rule=enforceRobustSUnc)
        mesg(settings['logFile'], '\tDone with steady state constraints\n')
    except:
        mesg(settings['logFile'], '\tError: Failed to define steady state constraints\n')
        return sol, stat

    try:
        # checks if there is any enzymes
        if min(Kspr.shape[0], Kspr.shape[1]) > 0:
            # Uncertain Kcat constraints
            if KUncFlag:
                # Flux enzyme relationships
                def ensureEnzymeBndKUnc(yst, i):
                    rowInd = numpy.nonzero(nonZeroIndK[0] == i)[0]
                    return sum(Kspr[i, nonZeroIndK[1][k]] * yst.v[nonZeroIndK[1][k]] for k in rowInd) \
                           - yst.e[P[i] - 1] + yst.Kw[i] <= bK[i]  # zero indexing, but permutation is 1 indexed

                yst.enzymeBndKUnc = Constraint(yst.enzK, rule=ensureEnzymeBndKUnc)

                # Equate Rx to u
                def equateRobustKUnc(yst, i):
                    rowInd = numpy.nonzero(nonZeroIndRK[0] == i)[0]
                    return sum(RKspr[i, nonZeroIndRK[1][k]] * yst.v[nonZeroIndRK[1][k]] for k in rowInd) - yst.Ku[
                        i] == 0

                yst.rbstBalKUnc = Constraint(yst.RKxid, rule=equateRobustKUnc)

                # Add the SOCP constraint
                def enforceRobustKUnc(yst, i):
                    return sum(yst.Ku[k] * yst.Ku[k] for k in range(RKindx[i], RKindx[i + 1])) - yst.Kw[i] * yst.Kw[
                        i] <= 0

                yst.socpConstKUnc = Constraint(yst.enzK, rule=enforceRobustKUnc)
            # Certain Kcat constraints
            else:
                def ensureEnzymeBnd(yst, i):
                    rowInd = numpy.nonzero(nonZeroIndK[0] == i)[0]
                    return sum(Kspr[i, nonZeroIndK[1][k]] * yst.v[nonZeroIndK[1][k]] for k in rowInd) \
                           - yst.e[P[i] - 1] <= bK[i]  # zero indexing, but permutation is 1 indexed

                yst.enzymeBnd = Constraint(yst.enzK, rule=ensureEnzymeBnd)
            mesg(settings['logFile'], '\tDone with enzyme-flux constraints\n')
    except:
        mesg(settings['logFile'], '\tError: Failed to define enzyme-flux constraints\n')
        return sol, stat

    try:
        # Pool constraint - bound enzyme mass
        def ensureEnzymeMassBnd(yst):
            # return (0,sum(MW[j] * yst.e[j] for j in yst.enzK),totMass)
            if all(MW[j] == 0 for j in range(len(MW))):
                return Constraint.Skip
            else:
                return sum(MW[j] * yst.e[j] for j in yst.enzK) <= totMass

        yst.enzymeMassBnd = Constraint(rule=ensureEnzymeMassBnd)
        mesg(settings['logFile'], '\tDone with the pool constraint\n')
    except:
        mesg(settings['logFile'], '\tError: Failed to define the pool constraints\n')
        return sol, stat

    try:
        mesg(settings['logFile'], '\tCreating model instance\n')
        instance = yst.create_instance()  # create a concrete instance of the model
        mesg(settings['logFile'], '\tDone creating model instance\n')
    except:
        mesg(settings['logFile'], '\tError: failed to create model instance\n')
        return sol, stat

    try:
        Opt = SolverFactory(settings['solver'])  # define our solver
        Opt.options["Crossover"] = settings['Crossover']  # A couple of options for Gurobi
        Opt.options["Method"] = settings['LPmethod']  # 0 = primal simplex, 1 = dual simplex, 2 = interior point
        Opt.options["NumericFocus"] = settings['NumericFocus']  # Values between 0 and 3 (default 0) with increasing accuracy
        Opt.options["BarQCPConvTol"] = settings['QPtol']  # conv tol for SOCP
        Opt.options["BarConvTol"] = settings['LPtol']  # conv tol for LP
        Opt.options['BarCorrectors'] = settings['BarrierCorrections']
        Opt.options['Presolve'] = settings['Presolve']  # -1 - auto, 0 - off, 1 - moderate, 2 - aggressive
        Opt.options['BarHomogeneous'] = settings['BarHomogeneous']  # alternative IP method that may better differentiate infeasibility and unboundedness
        Opt.options['Aggregate'] = settings['Aggregate']  # presolve parameter, might cause numerical troubles when on (=1, default value)
        start = time.time()
        Soln = Opt.solve(instance, tee=settings['solveOutput'])  # Solve the problem
        end = time.time()
        #instance.solutions.load_from(Soln)  # load the solution into the model TODO: verify that this is superfluous
        mesg(settings['logFile'], '\tDone solving the instance\n')
    except:
        mesg(settings['logFile'], '\tError: Failed to solve the instance\n')
        return sol, stat

    #
    # Check to make sure we solved the problem
    #
    if Soln.solver.termination_condition == TerminationCondition.optimal:
        # print "Optimal growth rate = "+str(value(instance.grwthRate))
        stat = 0
        mesg(settings['logFile'], '\tInstance was solved optimally\n')

        # only save enzyme solution if model has enzymes
        optE = []
        if min(Kspr.shape[0], Kspr.shape[1]) > 0:
            optE = [value(instance.e[i]) for i in range(Kspr.shape[0])]

        # save important values to solution as a dictionary
        sol = {'grwthRate': value(instance.grwthRate),
               'optFlux': [value(instance.v[i]) for i in range(SsprCer.shape[1])],
               'optE': optE,
               'SsprCer': SsprCer,
               'SsprUnc': SsprUnc,
               'cerRowIndS': cerRowInd,
               'uncerRowIndS': uncerRowInd,
               'cS': cS,
               'MW': MW,
               'lbS': lbS,
               'ubS': ubS,
               'solverTime': end-start}
        # save instance if told to
        if settings['retrieveInstance']:
            sol['instance'] = instance
        # we're done
        return sol, stat
    # elif Soln.status == SolverStatus.feasible: #TODO: is it possible to check for suboptimality?
    #     stat = 2  # stat = 2 is flag for suboptimality in uncModel()
    #     mesg(settings['logFile'], '\tInstance was solved suboptimally\n\n')
    #
    #     # only save enzyme solution if model has enzymes
    #     optE = []
    #     if min(Kspr.shape[0], Kspr.shape[1]) > 0:
    #         optE = [value(instance.e[i]) for i in range(Kspr.shape[0])]
    #
    #     # save important values to solution as a dictionary
    #     sol = {'grwthRate': value(instance.grwthRate),
    #            'optFlux': [value(instance.v[i]) for i in range(SsprCer.shape[1])],
    #            'optE': optE,
    #            'SsprCer': SsprCer,
    #            'SsprUnc': SsprUnc,
    #            'cerRowIndS': cerRowInd,
    #            'uncerRowIndS': uncerRowInd,
    #            'cS': cS}
    #     # save instance if told to
    #     if settings['retrieveInstance']:
    #         sol['instance'] = instance
    #     # we're done
    #     return sol, stat
    else:
        mesg(settings['logFile'], "Error: terminated with " + str(Soln.solver.termination_condition) + '\n\n')
        return sol, stat


