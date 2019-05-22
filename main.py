from input import Input
from methods import Methods
import pandas                                                                   
from scipy.stats import multivariate_normal                                     
import numpy as np                                                              

class Main(Input, Methods):
  def __init__(self):
    super().__init__()
    self.percent = []
    self.lrCrit = []
    self.lrValues = []
    self.lrPvalues = []
    self.lrMvalues = []
    self.lrSignif = []
    self.miCrit = []
    self.miValues = []
    self.miPvalues = []
    self.miMvalues = []
    self.miSignif = []
    self.null()
    self.analyze()
    self.output()

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
    lrLoc = []
    miLoc = []
    zs = self.combo[self.combo.columns[self.combo.columns.to_series().str.contains('Z_')]]
    zs = np.array(zs)
    self.lrValues = self.lr_real(zs)
    self.lrPvalues =  self.pvalues(self.lrCrit, self.lrValues, self.percent)
    self.lrSignif = np.where(self.lrPvalues<= 5e-8)
    self.lrSignif = np.squeeze(self.lrSignif, axis = 1)
    loc = self.combo[['CHR','BP']]
    #self.combo[self.combo.columns[self.combo.columns.to_series().str.contains("CHR|BP")]]
    for i in range(1,23):
      print('select: ' + str(i))
      select = loc[loc['CHR'] == i]
      print(select)
      select = loc['BP']
      print(select)
      select = np.array(select)
      print(select)
      lrLoc = select[self.lrSignif]
      print(lrLoc)
    exit()
    #print(loc)
    #if self.alpha is -1.0:
    #loc = self.combo[self.combo.columns[self.combo.columns.to_series().str.contains("CHR|BP")]]
    #print('lrLoc')
    #del loc
    self.lrMvalues = self.mvalues(zs[self.lrSignif,:], lrLoc)
    if self.migwas:
      self.miValues = self.mi_real(zs)
      self.miPvalues = self.pvalues(self.miCrit, self.miValues, self.percent)
      self.miSignif = np.where(self.miPvalues<= 5e-8)
      self.miSignif = np.squeeze(self.miSignif, axis = 1)
      #if self.alpha is -1.0:
      #loc = np.array(loc)
      print('miLoc')
      miLoc = loc[self.miSignif,:]
      #del loc
      self.miMvalues = self.mvalues(zs[self.miSignif,:], miLoc)

  #writes final output
  #TODO: needs to be updated to work, combo not actually created yet
  def output(self):
    if self.migwas:
      self.combo['MIGWAS_Score'] = self.miValues
      self.combo['MIGWAS_Pvalue'] = self.miPvalues
      for i in range(len(self.mean)):
        x = np.zeros(int(self.count))
        x[self.miSignif] = self.miMvalues[i]
        self.combo['MIGWAS_Mvalue_'+self.traits[i]] = x
    self.combo['PAT_Score'] = self.lrValues
    self.combo['PAT_Pvalue'] = self.lrPvalues
    for i in range(len(self.mean)):
      x = np.zeros(int(self.count))
      x[self.lrSignif] = self.lrMvalues[i]
      self.combo['PAT_Mvalue_'+self.traits[i]] = x
    self.combo.to_csv(self.out, sep=' ', index=False)
        
if __name__=="__main__":
  Main()
