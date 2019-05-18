from scipy.stats import multivariate_normal
import numpy as np

#Main class is responsible for input/output and for calling the two classes Migwas and Likelihood_Ratio which actually do the computation work.
class Pat():
  def __init__(self, mean, sigmaE, sigmaG):
    self.mean = mean
    self.sigmaE = sigmaE
    self.sigmaG = sigmaG

  #takes in the null sims, sampling weights and percents and returns critical values
  def process_null(self, sims, weigh, percent):
    null = multivariate_normal.logpdf(sims, self.mean, self.sigmaE)           
    alt = multivariate_normal.logpdf(sims, self.mean, np.add(self.sigmaE,self.sigmaG))             
    lrcrit = np.asarray(np.subtract(alt, null))                                             
    lrorder = lrcrit.argsort()
    lrcrit = np.stack((lrcrit, weigh))
    lrcrit = lrcrit[:,lrorder]
    lrcrit[1,] = np.cumsum(lrcrit[1,::-1])[::-1]
    lrcrit[1,] = np.divide(lrcrit[1,],total)
    lrkeep = lrcrit[1,].size - np.searchsorted(lrcrit[1,::-1], percent, side='left')
    lrcrit = lrcrit[0,]
    lrcrit = lrcrit[lrkeep]
    return lrcrit

  #pvalues maps the real data to critical values and returns the p-value
  def pvalues(self, lrcrit, lrvalues, percent):
    lrkeep = np.searchsorted(lrcrit, lrvalues, side='left')
    lrkeep = [int(i-1) if i > 0 else int(i) for i in lrkeep]
    lrpvalues = percent[lrkeep]
    return lrpvalues

if __name__=="__main__":
  Pat()
