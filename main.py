from inout import IO
from methods import Methods
import pandas                                                                   
from scipy.stats import multivariate_normal                                     
import numpy as np                                                              

class Main(IO, Methods):
  def __init__(self):
    super().__init__()
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
    sims = np.random.multivariate_normal(self.mean, self.impSigmaE, int(self.num))
    #this gets the lr of the sims, want null_cov/sim_cov to be used as a weight
    null = multivariate_normal.pdf(sims, self.mean, self.impSigmaE)
    alt = multivariate_normal.pdf(sims, self.mean, self.sigmaE)                
    weigh = np.asarray(alt/null)                                               
    self.lrCrit = self.lr_null(sims, weigh, self.percent)    
    
    if self.migwas:
      self.miCrit = self.mi_null(sims, weigh, self.percent)
        
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
    #self.combo = pandas.read_table(self.out)
    zs = self.combo[self.combo.columns[self.combo.columns.to_series().str.contains('Z_')]]
    zs = np.array(zs)
    print('zs: ' + str(zs))
    self.lrValues = self.lr_real(zs)
    self.lrPvalues =  self.pvalues(self.lrCrit, self.lrValues, self.percent)
    #self.lrMvalues = Mvalues.method(zs[self.signif,:])
    if self.migwas:
      self.miValues = self.mi_real(zs)
      self.miPvalues = self.pvalues(self.miCrit, self.miValues, self.percent)
      #self.miMvalues = Mvalues.method(zs[self.signif,:])
        
if __name__=="__main__":
  Main()
