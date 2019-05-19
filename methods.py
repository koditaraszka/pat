from scipy.stats import multivariate_normal
import math
import numpy as np

#This class contains the pat and migwas methods
class Methods():
  def __init__(self):
    self.mean = []
    self.sigmaE = []
    self.sigmaG = []

  #inheritance is weird
  def set_input(self, mean, sigmaE, sigmaG):
    self.mean = mean
    self.sigmaE = sigmaE
    self.sigmaG = sigmaG

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

  def mvalues(self, data):

#if __name__=="__main__":
#  Methods()
