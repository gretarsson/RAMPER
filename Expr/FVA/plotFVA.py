from matplotlib import pyplot as plt
import pickle


# load the files
with open('./Results/fvaS0.01-0.23', 'rb') as f:
    fvaCer = pickle.load(f)

with open('./Results/fvaS0.01-0.23+K', 'rb') as f:
    fvaUnc = pickle.load(f)

# initialize
yCerMax = []
yCerMin = []
yUncMax = []
yUncMin = []
yCerInt = []
yUncInt = []

# counting number of points that are skipped due to non-optimal termination conditions
nMissPnts = 0

# fill lists to plot
for i in range(len(fvaCer)-1):
    if fvaCer[i][0] != -1 and fvaUnc[i][0] != -1 and fvaCer[i][1] != -1 and fvaUnc[i][1] != -1:
        yCerMax.append(fvaCer[i][1])
        yCerInt.append(fvaCer[i][1] - fvaCer[i][0])
        # intervals
        yUncInt.append(fvaUnc[i][1] - fvaUnc[i][0])
        # max
        yUncMax.append(fvaUnc[i][1])
        yCerMin.append(fvaCer[i][0])
        yUncMin.append(fvaUnc[i][0])
    else:
        nMissPnts += 1

# print nunmber of missing points to screen
print 'Number of fluxes skipped due to accrued non-optimal termination conditions: '+str(nMissPnts)

# difference in intervals
diffInt = [yUncInt[i] - yCerInt[i] for i in range(len(yCerInt))]
#plt.scatter(range(len(yCerMax)), [yCerMin[i] - yUncMin[i] for i in range(len(yCerMax))], s=6)
plt.ylabel('Interval length difference')
plt.xlabel('Flux index')
plt.scatter(range(len(yCerInt)), diffInt, s=6)
#plt.scatter(range(len(yCerInt)), yUncInt, s=6)
plt.legend()
plt.show()