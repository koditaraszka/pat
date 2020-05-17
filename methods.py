from sklearn.cluster import KMeans
from scipy.stats import multivariate_normal, norm
import random
import math
import numpy as np
import itertools

#This class contains the pat and migwas methods
class Methods():

  #pvalues maps the real data to critical values and returns the p-value
  def pvalues(self, crit, values, percent):
    keep = np.searchsorted(crit, values, side='left')
    keep = [int(i-1) if i > 0 else int(i) for i in keep]
    return percent[keep]

  #helper for processing critical values
  def process_crit(self, crit, weigh, percent):
    total = np.sum(weigh)
    order = crit.argsort()
    crit = np.stack((crit, weigh))
    crit = crit[:,order]
    crit[1,] = np.cumsum(crit[1,::-1])[::-1]
    crit[1,] = np.divide(crit[1,],total)
    keep = crit[1,].size - np.searchsorted(crit[1,::-1], percent, side='left')
    crit = crit[0,]
    crit = crit[keep]
    return crit

  #helper for log of sums
  @staticmethod
  def sumlog(x, y):
    combine = np.stack((x,y), axis=1)
    high = np.amax(combine, axis=1)
    high = high.astype(dtype=np.float128)
    low = np.amin(combine, axis=1)
    low = low.astype(dtype=np.float128)
    diff = np.subtract(high,low)
    total = low + np.log1p(np.exp(diff))
    return total
  
  #interpretation
  def mvalues(self, zvals, loc, chrm, bp):
    k = zvals.shape[1]
    include = [list(i) for i in list(itertools.product([1.0, 0.0], repeat=k))]
    #alphaData = self.prune(loc, chrm, bp)
    #self.set_alpha(alphaData)
    print(self.alpha)
    sigmaG = (self.polygenic/self.alpha)*self.sigmaG
    mat = self.sigmaE + sigmaG
    value = multivariate_normal.logpdf(zvals, self.mean, mat)
    all_configs = np.copy(value)
    mvalues = [np.copy(value) for i in range(k)]
    for z in range(1,len(include)):
      alt = np.copy(mat)
      loc = include[z]
      for i in range(0,k):
        for j in range(0,i+1):
          if loc[i] == 0 or loc[j] == 0:
            alt[i,j] = self.sigmaE[i,j]
            alt[j,i] = self.sigmaE[j,i]
      value = multivariate_normal.logpdf(zvals, self.mean, alt)
      all_configs = Methods.sumlog(all_configs, value)
      for m in range(0,k):
        if loc[m] == 1.0:
          mvalues[m] = Methods.sumlog(mvalues[m],value)
    mvalues = np.exp(np.subtract(mvalues,all_configs))
    return mvalues

  #set alpha based on "independent" SNPs
  def prune(self, loc, chrm, bp):
    if self.sims:
      data=np.delete(loc,[chrm,bp], axis=1)
      print(data.shape)
      return data
    select = []
    for i in range(1,23):
      pruningChrm = loc[np.squeeze(np.where(loc[:,chrm] == i), axis = 1)]
      maxBP = np.amax(pruningChrm, axis = 0)[bp]
      top = (maxBP+100000) - (maxBP % 100000)
      groups = [float(i) for i in range(1, int(top), 100000)]
      for j in range(0, len(groups)-1):
        lower = groups[j]
        upper = groups[j+1]-1.0
        pruningBP = pruningChrm[np.squeeze(np.where((pruningChrm[:,bp] >= lower) & (pruningChrm[:,bp] <= upper)), axis = 1)]
        if pruningBP.ndim == 1:
          select.append(np.delete(pruningBP,[chrm,bp]))
        else:
          x=pruningBP.shape[0]
          if x > 0:
            keep = random.randint(0,(x-1))
            select.append(np.delete(pruningBP[keep,:],[chrm,bp]))
      lower = groups[-1]
      upper = maxBP
      pruningBP = pruningChrm[np.squeeze(np.where((pruningChrm[:,bp] >= lower) & (pruningChrm[:,bp] <= upper)), axis = 1)]
      if pruningBP.ndim == 1:
        select.append(np.delete(pruningBP,[chrm,bp]))
      else:
        x=pruningBP.shape[0]
        if x > 0:
          keep = random.randint(0,(x-1))
          select.append(np.delete(pruningBP[keep,:],[chrm,bp]))
    select = np.array(select)
    return select

  #set alpha
  def set_alpha(self, data):
    maximum = -np.inf
    for x in range(int(self.count),100,-25):
      if x == 0:
        x = 1
      x = float(x)
      mult = self.polygenic/x
      covar = self.sigmaE + mult*self.sigmaG
      pdf = multivariate_normal.logpdf(data, self.mean, covar)
      total = np.sum(pdf)
      if total > maximum:
        maximum = total
        self.alpha = x

  #returns critical values for likelihood ratio method
  def lr_null(self, sims, weigh, percent):
    null = multivariate_normal.logpdf(sims, self.mean, self.sigmaE)
    alt = multivariate_normal.logpdf(sims, self.mean, np.add(self.sigmaE,self.sigmaG))
    lrcrit = np.asarray(np.subtract(alt, null))
    return self.process_crit(lrcrit, weigh, percent)

  #returns lr results on real data
  def lr_real(self, data):
    null = multivariate_normal.logpdf(data, self.mean, self.sigmaE)
    alt = multivariate_normal.logpdf(data, self.mean, np.add(self.sigmaE,self.sigmaG))
    return np.asarray(np.subtract(alt, null))

  #return critical values for migwas method
  def mi_null(self, sims, weigh, percent):
    micrit = np.amax(abs(sims), axis=1)
    return self.process_crit(micrit, weigh, percent)

  #returns migwas results on real data
  def mi_real(self, data):
    return np.asarray(np.amax(abs(data), axis=1))
