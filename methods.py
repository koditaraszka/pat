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

  @staticmethod
  def sumlog(x, y):
    combine = np.stack((x,y), axis=1)
    high = np.amax(combine, axis=1)
    low = np.amin(combine, axis=1)
    diff = np.subtract(high,low)
    total = low + np.log1p(np.exp(diff))
    return total
  
  def mvalues(self, data, loc, chrm, bp):
    k = data.shape[1]
    include = [list(i) for i in list(itertools.product([1.0, 0.0], repeat=k))]
    groups = [[] for i in range(k)]
    mvalues = [[] for i in range(k)]
    configs = []
    #print('set alpha')
    print('alpha: ' + str(self.alpha))
    alphaData = self.prune(loc, chrm, bp)
    self.set_alpha(alphaData)
    for z in range(0,len(include)):
      alt = np.ones((k,k))
      loc = include[z]
      print('Loc: ' + str(loc))
      for i in range(0,k):
        for j in range(0,i+1):
          if loc[i] == 0 or loc[j] == 0:
            alt[i,j] = self.sigmaE[i,j]
            alt[j,i] = self.sigmaE[j,i]
          elif i == j:
            alt[i,i] = self.sigmaE[i,i] + self.alpha[i]*self.sigmaG[i,i]
          else:
            corrG = self.sigmaG[i,j]/(math.sqrt(self.sigmaG[i,i])*math.sqrt(self.sigmaG[j,j]))
            g = corrG*math.sqrt(self.alpha[i]*self.sigmaG[i,i])*math.sqrt(self.alpha[j]*self.sigmaG[j,j])
            alt[i,j] = self.sigmaE[i,j] + g
            alt[j,i] = self.sigmaE[j,i] + g
      configs.append(alt)
      for l in range(0,k):
        if loc[l] == 1.0:
          groups[l].append(z)
    all_configs = []
    for i in range(0,len(configs)):
      value = multivariate_normal.logpdf(data, self.mean, configs[i])
      if len(all_configs) is 0:
        all_configs = value
      else:
        all_configs = Methods.sumlog(all_configs, value)
      for j in range(0,len(groups)):
        if i in groups[j]:
          if len(mvalues[j]) is 0:
            mvalues[j] = value
          else:
            mvalues[j] = Methods.sumlog(mvalues[j], value)
    mvalues = np.exp(np.subtract(mvalues,all_configs))
    print('mvalues: ' + str(mvalues))
    return mvalues

  def prune(self, loc, chrm, bp):
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
    print(select)
    print('select: ' + str(select.shape))
    print('loc: ' + str(loc.shape))
    return select

  def set_alpha(self, data):
    k = data.shape[1]
    maximum = [-np.inf for i in range(k)]
    for i in range(100,100000,10):
      covar = self.sigmaE + (float(i)/100.0)*self.sigmaG
      for j in range(k):
        pdf = norm.logpdf(data[:,j], self.mean[j], covar[j,j])
        total = np.sum(pdf)
        if total > maximum[j]:
          maximum[j] = total
          self.alpha[j] = float(i)/100.0
    print('alpha: ' + str(self.alpha))

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
