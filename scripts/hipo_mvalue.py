import numpy as np
import pandas
import itertools
from scipy.stats import multivariate_normal, norm


def sumlogs(x,y):
  combine = np.stack((x,y), axis=1)
  high = np.amax(combine, axis=1)
  high = high.astype(dtype=np.float128)
  low = np.amin(combine, axis=1)
  low = low.astype(dtype=np.float128)
  diff = np.subtract(high,low)
  total = low + np.log1p(np.exp(diff))
  return total

combo = pandas.read_csv("real_analysis.txt", sep ='\t')
combo['chrm'] = 1
loc = combo[combo.columns[combo.columns.to_series().str.contains("chr|rsid|z.")]]
del loc["z.HIPOD1"]
del loc["z.HIPOD2"]
del loc["z.HIPOD3"]
del loc["z.HIPOD4"]

loc = np.array(loc)
signif = combo.loc[(combo["pval.HIPOD1"] <= 5e-8) | (combo["pval.HIPOD2"] <= 5e-8) | (combo["pval.HIPOD3"] <= 5e-8) | (combo["pval.HIPOD4"] <= 5e-8)].index.values

zs = combo[combo.columns[combo.columns.to_series().str.contains('z.')]]
del zs["z.HIPOD1"]
del zs["z.HIPOD2"]
del zs["z.HIPOD3"]
del zs["z.HIPOD4"]

zs = np.array(zs)

k = 4
include = [list(i) for i in list(itertools.product([1.0, 0.0], repeat=k))]

alphaData = np.delete(loc,[0,5], axis = 1)
alphaData = alphaData[signif,:]
print(alphaData.shape)
count = 7025734

sigmaG = np.array([[ 0.06248276,0.01365356,-0.01424399,0.00799307],
    [ 0.01365356,0.03207242,-0.0075916,0.02182045],
    [-0.01424399,-0.0075916,0.12072998,-0.00964204],
    [ 0.00799307,0.02182045,-0.00964204,0.03377293]])

sigmaE = np.array( [[ 1.00000000e+00,2.21634186e-01,-6.59629045e-02,1.45811322e-01],
    [ 2.21634186e-01,1.00000000e+00,-9.71533404e-04,6.77997010e-01],
    [-6.59629045e-02,-9.71533404e-04,1.00000000e+00,-4.56618686e-02],
    [ 1.45811322e-01,6.77997010e-01,-4.56618686e-02,1.00000000e+00]])


maximum = -np.inf
alpha = -1

for x in range(int(count),100,-25):
  if x == 0:
    x = 1
  x = float(x)
  mult = count/x
  covar = sigmaE + mult*sigmaG
  pdf = multivariate_normal.logpdf(alphaData,np.array([0.0,0.0,0.0,0.0]), covar)
  total = np.sum(pdf)
  if total > maximum:
    maximum = total
    alpha = x

print(alpha)
sigmaG = (count/alpha)*sigmaG
mat = sigmaE + sigmaG
value = multivariate_normal.logpdf(zs, np.array([0.0,0.0,0.0,0.0]), mat)
all_configs = np.copy(value)
mvalues = [np.copy(value) for i in range(k)]
for z in range(1,len(include)):
  alt = np.copy(mat)
  loc = include[z]
  for i in range(0,k):
    for j in range(0,i+1):
      if loc[i] == 0 or loc[j] == 0:
        alt[i,j] = sigmaE[i,j]
        alt[j,i] = sigmaE[j,i]
  value = multivariate_normal.logpdf(zs, np.array([0.0,0.0,0.0,0.0]), alt)
  all_configs = sumlogs(all_configs, value)
  for m in range(0,k):
    if loc[m] == 1.0:
      mvalues[m] = sumlogs(mvalues[m],value)
mvalues = np.exp(np.subtract(mvalues,all_configs))

for i in range(k):
  x = np.zeros(int(count))
  x = mvalues[i]
  combo['HIPO_Mvalue_'+str(i)] = x
    
combo.to_csv("real_mvalue.txt", sep='\t', index=False)
