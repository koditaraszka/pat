from inout import IO
from migwas import Migwas
from pat import Pat
import pandas                                                                   
from scipy.stats import multivariate_normal                                     
import numpy as np                                                              

class Main():
  inOut = IO()
  if inOut.migwas:
    mi = Migwas()
  lr = Pat(inOut.mean, inOut.sigmaE, inOut.sigmaG)
  def __init__(self):
    self.percent = []
    self.count = 0
    self.lrCrit = []
    self.lrValues = []
    self.lrPvalues = []
    self.lrMvalues = []
    
    self.miCrit = []
    self.miValues = []
    self.miPvalues = []
    self.miMvalues = []
    
    self.signif = 0
    self.null()
    self.analyze()

  #null actually simulates the data
  def null(self):
    self.set_percent()
    sims = np.random.multivariate_normal(inOut.mean, inOut.impSigmaE, int(inOut.num))
    #this gets the lr of the sims, want null_cov/sim_cov to be used as a weight
    null = multivariate_normal.pdf(sims, inOut.mean, inOut.impSigmaE)
    alt = multivariate_normal.pdf(sims, inOut.mean, inOut.sigmaE)                
    weigh = np.asarray(alt/null)                                               
    total = np.sum(weigh)
    self.lrCrit = lr.process_null(sims, weigh, self.percent)    
    
    if inOut.migwas:
      mi = Migwas()
      self.miCrit = mi.process_null(sims, weigh, self.percent)
        
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
    
    #if smaller map 1.0 & if larger map to 1e-8
    percent = [1.0] + x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8
    self.percent = np.asarray(percent)
        
  #analyze begins the actual analysis of the real input
  #TODO Needs to be generalized to different input.
  def analyze(self):
    inOut.combo = pandas.read_table(inOut.out)
    zs = inOut.combo[inOut.combo.columns[inOut.combo.columns.to_series().str.contains('Z_')]]
    zs = np.array(zs)
    self.lrValues = lr.process_real(zs)
    self.lrPvalues =  lr.pvalues(self.lrCrit, self.lrValues, self.percent)
    #self.lrMvalues = Mvalues.method(zs[self.signif,:])
    if inOut.migwas:
      self.miValues = mi.process_real(zs)
      self.miPvalues = mi.pvalues(self.miCrit, self.miValues, self.percent)
      #self.miMvalues = Mvalues.method(zs[self.signif,:])
        
if __name__=="__main__":                                                        
  Main()  
