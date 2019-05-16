class Main():
  def __init__(self):
    init = IO()  
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
    self.null(init)
    self.set_percent()
    self.analyze()

  #null actually simulates the data
  def null(self,init):
    sims = np.random.multivariate_normal(init.mean, init.impSigmaE, int(self.num))
    #this gets the lr of the sims, want null_cov/sim_cov to be used as a weight
    w = Likelihood_Ratio(sims, init.mean, init.impSigmaE, init.sigmaE, [], [], log=False)
    weigh = w.get_values()
    total = np.sum(weigh)

    #need to break things into parts because had memory issues previously

    #this runs the likelihood ratio method on simulations
    lr = Likelihood_Ratio(sims, init.mean, init.sigmaE, init.SigmaG, [], [])
    #this runs migwas on the simulations
    if init.migwas:
      mi = Migwas(sims)
    del sims
    
    #these are the unsorted critical values
    self.lrcrit = lr.get_values()
    lrorder = self.lrcrit.argsort()
    #maps the weighing of the critical value to the critical value
    self.lrcrit = np.stack((self.lrcrit, weigh))
    #sorts the critical values by the original ordering
    self.lrcrit = self.lrcrit[:,lrorder]
    del lrorder
    
    if init.migwas:
      self.micrit = mi.get_values()
      miorder = self.micrit.argsort()
      self.micrit = np.stack((self.micrit, weigh))
      self.micrit = self.micrit[:,miorder]
      del miorder
    del weigh
    
    
    #this now sums all weights up to this location and then divides by total weight
    self.lrcrit[1,] = np.cumsum(self.lrcrit[1,::-1])[::-1]
    self.lrcrit[1,] = np.divide(self.lrcrit[1,],total)
    if init.migwas:
      self.micrit[1,] = np.cumsum(self.micrit[1,::-1])[::-1]
      self.micrit[1,] = np.divide(self.micrit[1,],total)
        
  #set_percent starts matching critical values to a p-value (1-alpha)
  def set_percent(self):
    w = list(range(9999,999,-1))
    w = [round(float(i)/(1e4),4) for i in w]
    #1e-1
    x1 = [round(float(i)*1e0,4) for i in w]
    #1e-2
    x2 = [round(float(i)*1e-1,5) for i in w]
    #1e-3
    x3 = [round(float(i)*1e-2,6) for i in w]
    #1e-4
    x4 = [round(float(i)*1e-3,7) for i in w]
    #1e-5
    x5 = [round(float(i)*1e-4,8) for i in w]
    #1e-6
    x6 = [round(float(i)*1e-5,9) for i in w]
    #1e-7
    x7 = [round(float(i)*1e-6,10) for i in w]
    #1e-8
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

    if init.migwas:
      mikeep = self.micrit[1,].size - np.searchsorted(self.micrit[1,::-1],self.percent,side='left')
      self.micrit = self.micrit[0,]
      self.micrit = self.micrit[mikeep]
        
  #analyze begins the actual analysis of the "real" input data
  #TODO Needs to be generalized to different input.
  def analyze(self, init):
    init.combo = pandas.read_table(init.out)
    zs = init.combo[init.combo.columns[init.combo.columns.to_series().str.contains('Z_')]]
    zs = np.array(zs)
    m = Migwas(zs)
    self.mivalues = m.get_values()
    l = Likelihood_Ratio(zs, init.mean, init.sigmaE, init.SigmaG, init.cfigSigmaG, init.groups)
    self.lrvalues = l.get_values()
    self.pvalues()
    self.mvalues = l.set_mvalues(zs[self.signif,:])
        
        
  #pvalues maps the real data to the corresponding critical value and it's pvalue
  def pvalues(self,init):
    lrkeep = np.searchsorted(self.lrcrit,self.lrvalues, side='left')
    lrkeep = [int(i-1) if i > 0 else int(i) for i in lrkeep]
    self.lrpvalues = self.percent[lrkeep]
    if init.migwas:
      mikeep = np.searchsorted(self.micrit,self.mivalues, side='left')
      mikeep = [int(i-1) if i > 0 else int(i) for i in mikeep]
      self.mipvalues = self.percent[mikeep]
    self.signif = np.where(self.lrpvalues<= 5e-8)
    self.signif = np.squeeze(self.signif, axis = 1)
    
