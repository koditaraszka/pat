from in_out import IO
#import numpy as  np
#import Migwas
import numpy.linalg as linalg                                                   
import pandas                                                                   
from scipy.stats import multivariate_normal                                     
import numpy as np                                                              
import argparse                                                                 
import math                                                                     
import itertools

class Pat():
  inOut = IO()  
  def __init__(self):
    self.percent = []
    self.count = 0
    self.lrcrit = []
    self.lrvalues = []
    self.lrpvalues = []
    self.micrit = []
    self.mivalues = []
    self.mipvalues = []
    self.mvalues = []
    self.signif = 0
    self.null()
    self.set_percent()
    self.analyze()

  #null actually simulates the data
  def null(self):
    sims = np.random.multivariate_normal(inOut.mean, inOut.impSigmaE, int(inOut.num))
    #this gets the lr of the sims, want null_cov/sim_cov to be used as a weight
    null = multivariate_normal.pdf(sims, inOut.mean, inOut.impSigmaE)              
    alt = multivariate_normal.pdf(sims, inOut.mean, inOut.sigmaE)                
    weigh = np.asarray(alt/null)                                               
    total = np.sum(weigh)

    #need to break things into parts because had memory issues previously

    #this runs the likelihood ratio method on simulations
    self.lrcrit = self.pat(sims, inOut.mean, inOut.sigmaE, inOut.SigmaG, [], [])
    #this runs migwas on the simulations
    if inOut.migwas:
      mi = Migwas(sims)
    del sims
    
    #these are the unsorted critical values
    lrorder = self.lrcrit.argsort()
    #maps the weighing of the critical value to the critical value
    self.lrcrit = np.stack((self.lrcrit, weigh))
    #sorts the critical values by the original ordering
    self.lrcrit = self.lrcrit[:,lrorder]
    del lrorder
    
    if inOut.migwas:
      self.micrit = mi.get_values()
      miorder = self.micrit.argsort()
      self.micrit = np.stack((self.micrit, weigh))
      self.micrit = self.micrit[:,miorder]
      del miorder
    del weigh
    
    #this now sums all weights up to this location and then divides by total weight
    self.lrcrit[1,] = np.cumsum(self.lrcrit[1,::-1])[::-1]
    self.lrcrit[1,] = np.divide(self.lrcrit[1,],total)
    if inOut.migwas:
      self.micrit[1,] = np.cumsum(self.micrit[1,::-1])[::-1]
      self.micrit[1,] = np.divide(self.micrit[1,],total)
        
  #set_percent starts matching critical values to a p-value (1-alpha)
  def set_percent(self):
    w = list(range(9999,999,-1))
    w = [round(float(i)/(1e4),4) for i in w]
    x1 = [round(float(i)*1e0,4) for i in w]
    x2 = [round(float(i)*1e-1,5) for i in w]
    x3 = [round(float(i)*1e-2,6) for i in w]
    x4 = [round(float(i)*1e-3,7) for i in w]
    x5 = [round(float(i)*1e-4,8) for i in w]
    x6 = [round(float(i)*1e-5,9) for i in w]
    x7 = [round(float(i)*1e-6,10) for i in w]
    x8 = [round(float(i)*1e-7,11) for i in w]
    
    #values smaller than the smallest critical value will map to 1.0
    #values larger than our largest critical value will map to 1e-8
    percent = [1.0] + x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8
    self.percent = np.asarray(percent)
    self.map_values()
        
  #map_values take the sum(weigh[i:])/total that is >= to the percent of interest
  def map_values(self):
    lrkeep = self.lrcrit[1,].size - np.searchsorted(self.lrcrit[1,::-1],self.percent,side='left')
    #only keep the critical vaules
    self.lrcrit = self.lrcrit[0,]
    #only keep the critical values that mapped
    self.lrcrit = self.lrcrit[lrkeep]

    if inOut.migwas:
      mikeep = self.micrit[1,].size - np.searchsorted(self.micrit[1,::-1],self.percent,side='left')
      self.micrit = self.micrit[0,]
      self.micrit = self.micrit[mikeep]
        
  #analyze begins the actual analysis of the "real" input data
  #TODO Needs to be generalized to different input.
  def analyze(self):
    inOut.combo = pandas.read_table(inOut.out)
    zs = inOut.combo[inOut.combo.columns[inOut.combo.columns.to_series().str.contains('Z_')]]
    zs = np.array(zs)
    m = Migwas(zs)
    self.mivalues = m.get_values()
    self.lrvalues = self.pat(zs, inOut.mean, inOut.sigmaE, inOut.SigmaG, inOut.cfigSigmaG, inOut.groups)
    self.pvalues()
    #self.mvalues = l.set_mvalues(zs[self.signif,:])
        
  #pvalues maps the real data to the corresponding critical value and it's pvalue
  def pvalues(self):
    lrkeep = np.searchsorted(self.lrcrit,self.lrvalues, side='left')
    lrkeep = [int(i-1) if i > 0 else int(i) for i in lrkeep]
    self.lrpvalues = self.percent[lrkeep]
    if inOut.migwas:
      mikeep = np.searchsorted(self.micrit,self.mivalues, side='left')
      mikeep = [int(i-1) if i > 0 else int(i) for i in mikeep]
      self.mipvalues = self.percent[mikeep]
    self.signif = np.where(self.lrpvalues<= 5e-8)
    self.signif = np.squeeze(self.signif, axis = 1)

  #likelihood ratio for test statistic                                                
  def pat(self, sims, mean, null_cov, alt_cov, alt_configs, groups):
    null = multivariate_normal.logpdf(sims, mean, null_cov)         
    alt = multivariate_normal.logpdf(sims, mean, alt_cov)           
    ratio = alt - null                                                        
    return np.asarray(ratio)
                                                                               
if __name__=="__main__":                                                        
    Pat()  
