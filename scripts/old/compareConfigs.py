'''
Author Kodi Collins
Email: kodicollins@ucla.edu
This script performs GWAS by leveraging genetic correlation between traits.
Run python theMethod.py -h or see the readme for more information on how to use the method.  
There is some hardcoding in this method:
  1. the trait/environmental variance-covariance is assumed to be already standardized 
'''
import numpy.linalg as linalg
import bisect
import pandas
from scipy.stats import multivariate_normal
import numpy as np
import argparse
import math
import itertools

#Main class is responsible for input/output and for calling the two classes Migwas and Likelihood_Ratio which actually do the computation work.
class Main():
  def __init__(self):
    self.mean = []
    self.alt_cov = []
    self.null_cov = []
    self.sim_cov = []
    self.imp_var = 5.0
    self.percent = []
    self.count = 0
    self.size = []
    self.file = []
    self.gen = []
    self.env = []
    self.lrcrit = []
    self.lrvalues = []
    self.lrpvalues = []
    self.micrit = []
    self.mivalues = []
    self.mipvalues = []
    self.combo = []
    self.include = []
    self.alt_configs = []
    self.k = 0
    self.mvalues = []
    self.full_prob = []
    self.groups = []
    self.num = 1e6
    self.signif = 0
    self.outfile = 'gwas_results.txt'
    self.divide = 1.0
    self.define_parser()
    self.simulate()
    self.set_percent()
    self.analyze()
    self.output()

  #define_parser is the command line parser 
  def define_parser(self):
    parser = argparse.ArgumentParser(description='This program computes critical values, p-values and m-values for GWAS in multiple traits.')
    #required
    required = parser.add_argument_group('Required Arguments')
    required.add_argument('-e', '--envir', dest = 'envir', required = True, 
      help = 'text file where each line is as such trait names pairwise environmental correlation whitespace separated (eg. bmi height -0.088177 which is equivalent to height bmi -0.088177).')
    required.add_argument('-g', '--genetic', dest = 'genetic', required = True,
      help = 'text file where each line is the name of trait i and trait j, the genetic variance for trait i, genetic variance for trait j, and their pairwise genetic correlation (eg. bmi height 0.2468 0.453 -0.1589 which is equivalent to height bmi 0.453 0.2468 -0.1589)')
    required.add_argument('-n','--pop_sizes', dest = 'pop_sizes', required = True,
    help='text file where each line is the names of two traits and the sample size for trait 1, sample size for trait 2, and the overlapping individuals (eg. bmi height 336107 336474 336107 which is eqivalent to height bmi 336474 336107 336107')
    required.add_argument('-f', '--gwas_files', dest = 'gwas_files', required = True,
      help = 'text file where each line is the name of a trait and its corresponding GWAS summary statistics file (eg. bmi bmi.txt)')
    #optional
    optional = parser.add_argument_group('Optional Arguments')
    optional.add_argument('-s', '--sampling', dest = 'sampling', required = True,
      help = 'This takes in an number, which will be used as the variance in the the null simulations. The default is 5.')
    optional.add_argument('-x', '--num', dest = 'num', help = 'number of simulations to run. The default is 1e6.')
    optional.add_argument('-o', '--out', dest = 'out', help = 'text file where the output will be written. The default is gwas_results.txt')
    optional.add_argument('-d', '--divide', dest = 'divide', help = 'for mvalues')
    args = parser.parse_args()
    self.read_parser(args)

  #read_parser actually reads the command line arguments and is directly called from define_parser
  #while not exhaustive, the method does check for some invalid arguments
  #it checks that -envir, -genetic, -pop_sizes, -gwas_files are actually called
  #it checks that the input sizes are correct in relation to each other, this is not exhaustive
  #it checks that -model full (the default model) and the -prior flag are not used together
  def read_parser(self, args):
    #set the required arguments
    self.env = [line.split() for line in open(args.envir)]
    self.gen = [line.split() for line in open(args.genetic)]
    self.size = [line.split() for line in open(args.pop_sizes)]
    self.file = [line.split() for line in open(args.gwas_files)]
    self.k = len(self.file)
    self.groups = [[] for i in range(self.k)]
    self.imp_var = float(args.sampling)
    #set the other arguments if called
    if args.out is not None:
      self.outfile = args.out
    if args.divide is not None:
      self.divide = float(args.divide)
    if args.num is not None:
      self.num = int(args.num)
    if self.size is None or self.env is None or self.gen is None or self.file is None:
      raise ValueError("-envir, -genetic, -pop_sizes, and -gwas_files are required inputs, at least one is missing")
    if len(self.size) != len(self.env) and len(self.size) != len(self.gen) and len(self.size) != len(self.file):
      raise ValueError("Number of Pairwise Comparisons are not equal across input files")
    if (self.k*(self.k - 1)/2.0) != len(self.size):
      raise ValueError("k(k-1)/2 for the k traits from gwas_files does not equal the number of lines of pairwise comparisons in the -pop_size input file")

  #simulate sets the parameters of the models.
  def simulate(self):
    #2^k models for k traits
    self.include = [list(i) for i in list(itertools.product([1.0, 0.0],repeat=self.k))]
    #get the number of SNPs to change genetic variance to per snp
    self.set_persnp()
    e = np.ones((self.k, self.k))
    g = np.ones((self.k, self.k))
    p = np.ones((self.k, self.k))
    sigmaG = np.ones((self.k, self.k))
    sigmaE = np.ones((self.k, self.k))
    sim_v = np.ones((self.k, self.k))
    count = dict()
    num = 0
    for j in self.file:
      count[j[0]] = num
      num += 1
    
    #this puts the data into the list of lists of values.
    for i in range(0,len(self.env)):
      env_name1 = self.env[i][0]
      env_name2 = self.env[i][1]
      env_corr = float(self.env[i][2])
      if env_corr < -1 or env_corr > 1:
        raise ValueError("The environmenantal correlation value is either less than -1 or greater than 1, did you pass in a non-standardized environmental covariance?")
      one = count[env_name1]
      two = count[env_name2]
     
      #environmental covariance
      e[one][one] = 1.0
      e[two][two] = 1.0
      e[one][two] = env_corr
      e[two][one] = env_corr
      
      #simulated covariance
      sim_v[one][one] = self.imp_var
      sim_v[two][two] = self.imp_var
      sim_v[one][two] = env_corr*self.imp_var
      sim_v[two][one] = env_corr*self.imp_var

      #genetic covariance
      gen_name1 = self.gen[i][0]
      gen_name2 = self.gen[i][1]
      gen_var1 = float(self.gen[i][2])
      gen_var1 = gen_var1/float(self.count)
      gen_var2 = float(self.gen[i][3])
      gen_var2 = gen_var2/float(self.count)
      gen_corr = float(self.gen[i][4])
      gen_cov = gen_corr*math.sqrt(gen_var1)*math.sqrt(gen_var2)
      
      one = count[gen_name1]
      two = count[gen_name2]      
      g[one][one] = gen_var1
      g[two][two] = gen_var2
      g[one][two] = gen_cov
      g[two][one] = gen_cov

      pop_name1 = self.size[i][0]
      pop_name2 = self.size[i][1]
      pop_size1 = float(self.size[i][2])
      pop_size2 = float(self.size[i][3])
      pop_shared = float(self.size[i][4])
      pop_ratio = pop_shared/(math.sqrt(pop_size1)*math.sqrt(pop_size2))
      one = count[pop_name1]
      two = count[pop_name2]
      p[one][one] = pop_size1
      p[two][two] = pop_size2
      p[one][two] = pop_ratio
      p[two][one] = pop_ratio
    
    #This sets the sigmaG and sigmaE values
    #TODO go between 0 and math.ceil(self.k/2), but meh 
    for i in range(0,self.k):
      for j in range(0,self.k):
        if i == j:
          sigmaE[i][i] = e[i][i]
          sigmaG[i][i] = p[i][i]*g[i][i]
        else:
          sigmaE[i][j] = p[i][j]*e[i][j]
          sigmaE[j][i] = sigmaE[i][j]
          sigmaG[i][j] = math.sqrt(p[i][i])*math.sqrt(p[j][j])*g[i][j]
          sigmaG[j][i] = sigmaG[i][j]

    #null_cov is just the environmental covariance matrix
    self.null_cov = sigmaE
    #sim_cov is the actual "environmental covariance matrix" 
    self.sim_cov = sim_v
    
    #setting up alt config for mvalues and assigning index to group
    for z in range(0,len(self.include)):
      alt = np.zeros((self.k,self.k))
      all_alt = np.zeros((self.k,self.k))
      i = self.include[z]
      if 1.0 in i or 0.0 in i: #nonesense
        for j in range(0,self.k):
          for k in range(0,j+1):
            if i[j] == 0 or i[k] == 0:
              alt[j,k] = sigmaE[j,k]
              alt[k,j] = sigmaE[k,j]
            else:
              alt[k,j] = self.divide*sigmaG[k,j] + sigmaE[k,j]
              alt[j,k] = self.divide*sigmaG[j,k] + sigmaE[j,k]
            if 0.0 not in i:
              all_alt[k,j] = sigmaG[k,j] + sigmaE[k,j]
              all_alt[j,k] = sigmaG[j,k] + sigmaE[j,k]
        #print('i: ' + str(i))
        #print(alt)
        self.alt_configs.append(alt)
        if 0.0 not in i:
          self.alt_cov = all_alt
        #assign to groups
        for l in range(0,self.k):
          if i[l] == 1.0:
            self.groups[l].append(z)
    
    #print('all_alt')
    #print(self.alt_cov)
    #The mean is a vector of k 0s
    self.mean = np.zeros(self.k)
    print('self.count: ' + str(self.count))
    print('self.divide: ' + str(self.divide))
    self.null()


  #null actually simulates the data
  def null(self):
    sims = np.random.multivariate_normal(self.mean, self.sim_cov, int(self.num))
    #this gets the lr of the sims, want null_cov/sim_cov to be used as a weight
    w = Likelihood_Ratio(sims, self.mean, self.sim_cov, self.null_cov, [], [], log=False) 
    weigh = w.get_values()
    total = np.sum(weigh)
    
    #this runs the likelihood ratio method on simulations
    lr = Likelihood_Ratio(sims, self.mean, self.null_cov, self.alt_cov, [], [])
    #this runs migwas on the simulations
    mi = Migwas(sims)
    del sims
    
    #these are the unsorted critical values
    self.lrcrit = lr.get_values()
    self.micrit = mi.get_values()
    lrorder = self.lrcrit.argsort()
    miorder = self.micrit.argsort()
    
    #maps the weighing of the critical value to the critical value
    self.lrcrit = np.stack((self.lrcrit, weigh))
    self.micrit = np.stack((self.micrit, weigh))
    del weigh
    
    #sorts the critical values by the original ordering
    self.lrcrit = self.lrcrit[:,lrorder]
    self.micrit = self.micrit[:,miorder]
    del lrorder
    del miorder
    
    #this now sums all weights up to this location and then divides by total weight
    self.lrcrit[1,] = np.cumsum(self.lrcrit[1,::-1])[::-1]
    self.micrit[1,] = np.cumsum(self.micrit[1,::-1])[::-1]
    self.lrcrit[1,] = np.divide(self.lrcrit[1,],total)
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
    mikeep = self.micrit[1,].size - np.searchsorted(self.micrit[1,::-1],self.percent,side='left')
    
    #only keep the critical vaules
    self.lrcrit = self.lrcrit[0,]
    self.micrit = self.micrit[0,]
    #only keep the critical values that mapped
    self.lrcrit = self.lrcrit[lrkeep]
    self.micrit = self.micrit[mikeep]
    
  #analyze begins the actual analysis of the "real" input data
  #TODO this will need to be generalized to not different headers. 
  def analyze(self):
    traits = []
    files = []
    for i in self.file:
      traits.append(i[0])
      files.append(i[1])
    gwas = pandas.read_table(files[0])
    gwas = gwas.add_suffix('_'+traits[0])
    gwas = gwas.rename(index=str, columns = {'SNP_'+traits[0]: 'SNP', 'CHR_'+traits[0]: 'CHR', 'BP_'+traits[0]: 'BP', 'A1_'+traits[0]: 'A1', 'A2_'+traits[0]: 'A2'})
    for i in range(1,self.k):
      t = pandas.read_table(files[i])
      t = t.add_suffix('_'+traits[i])
      t = t.rename(index=str, columns = {'SNP_'+traits[i]: 'SNP', 'CHR_'+traits[i]: 'CHR', 'BP_'+traits[i]: 'BP', 'A1_'+traits[i]: 'A1', 'A2_'+traits[i]: 'A2'})
      gwas = gwas.merge(t, how = 'inner', on=['SNP','CHR','BP','A1','A2'])
    self.combo = gwas
    zs = gwas[gwas.columns[gwas.columns.to_series().str.contains('Z_')]]
    zs = np.array(zs)
    m = Migwas(zs)
    self.mivalues = m.get_values()
    l = Likelihood_Ratio(zs, self.mean, self.null_cov, self.alt_cov, self.alt_configs, self.groups)
    self.lrvalues = l.get_values()
    self.pvalues()
    print('number of zs: ' + str(zs.shape))
    #zs = np.transpose(zs)
    print('zs[0]: ' + str(zs[0]))
    print('zs[:,0]: ' + str(zs[:,0]))
    #zs = zs[:,self.signif]
    self.mvalues = l.set_mvalues(zs[self.signif,:])
    self.full_prob = l.get_fullprob(zs[self.signif,:])
  
  #set_persnp gets the per SNP heritability
  def set_persnp(self):
    traits = []
    files = []
    Z = []
    P = []
    for i in self.file:
      traits.append(i[0])
      files.append(i[1])
      Z.append('Z_'+i[0])
      P.append('P_'+i[0])
    gwas = pandas.read_table(files[0])
    gwas = gwas.add_suffix('_'+traits[0])
    gwas = gwas.rename(index=str, columns = {'SNP_'+traits[0]: 'SNP', 'CHR_'+traits[0]: 'CHR', 'BP_'+traits[0]: 'BP', 'A1_'+traits[0]: 'A1', 'A2_'+traits[0]: 'A2'})
    for i in range(1,self.k):
      t = pandas.read_table(files[i])
      t = t.add_suffix('_'+traits[i])
      t = t.rename(index=str, columns = {'SNP_'+traits[i]: 'SNP', 'CHR_'+traits[i]: 'CHR', 'BP_'+traits[i]: 'BP', 'A1_'+traits[i]: 'A1', 'A2_'+traits[i]: 'A2'})
      gwas = gwas.merge(t, how = 'inner', on=['SNP','CHR','BP','A1','A2'])
    self.count = float(gwas.shape[0])
    if self.divide != 1.0:
      self.divide = (self.count/self.divide)
    del gwas

  #pvalues maps the real data to the corresponding critical value and it's pvalue
  def pvalues(self):
    mikeep = np.searchsorted(self.micrit,self.mivalues, side='left')
    lrkeep = np.searchsorted(self.lrcrit,self.lrvalues, side='left')
    mikeep = [int(i-1) if i > 0 else int(i) for i in mikeep]
    lrkeep = [int(i-1) if i > 0 else int(i) for i in lrkeep]
    self.mipvalues = self.percent[mikeep]
    self.lrpvalues = self.percent[lrkeep]
    print(self.lrpvalues)
    print(' ')
    self.signif = np.where(self.lrpvalues<= 5e-8)
    self.signif = np.squeeze(self.signif, axis = 1)
    print(self.signif)
    print('self.signif: ' + str(self.signif.shape))

  #output now adds the MIGWAS and LR information to the a table of all of the traits
  def output(self):
    self.combo['MIGWAS_Score'] = self.mivalues
    self.combo['MIGWAS_Pvalue'] = self.mipvalues
    self.combo['LR_Score'] = self.lrvalues
    self.combo['LR_Pvalue'] = self.lrpvalues
    #for i in range(len(self.file)):
    x = np.zeros(int(self.count))
    x[self.signif] = self.mvalues
    self.combo['maxpvalue'] = x
    x = np.zeros(int(self.count))
    x[self.signif] = self.full_prob
    self.combo['fullpvalue'] = x
    self.combo.to_csv(self.outfile, sep=' ', index=False)

#Migwas is the class for computing the MIGWAS test statistic
class Migwas:
  #MIGWAS only needs the zscores
  def __init__(self,sims):
    self.results = []
    self.max_math(sims)

  #get_values returns the test statistic, but it does not sort them
  def get_values(self):
    return np.asarray(self.results)

  #max_math does the actual computation of test statistic
  def max_math(self, sims):
    self.results = np.amax(abs(sims), axis=1)

#Likelihood_Ratio is the class that computes the likelihood ratio test statistic
class Likelihood_Ratio():
  def __init__(self, sims, mean, null_cov, alt_cov, alt_configs, groups, log=True):
    print('sims shape: ' + str(sims.shape))
    self.mean = mean
    self.null_cov = null_cov
    self.alt_cov = alt_cov
    self.alt_configs = alt_configs
    self.groups = groups
    self.results = []
    self.mvalues = [[] for i in range(len(self.mean))]
    if log is False:
      self.ratio_math(sims)
    else:
      self.ratio_log_math(sims)
  
  #how to compute log mvalues
  @staticmethod                                                                 
  def sumlog(x, y):
    #print('x: ' + str(x))
    #print('y: ' + str(y))
    combine = np.stack((x,y), axis=1)
    #print('combine: ' + str(combine))
    #high = np.amax(combine, axis=1)
    #return high
    #print('high: ' + str(high))
    low = np.mean(combine, axis=1)
    return low
    #print('low: ' + str(low))
    #diff = np.subtract(high,low)
    #print('diff: ' + str(diff))
    #total = low + np.log1p(np.exp(diff))                 
    #print('total: ' + str(total))
    #return total

  #ratio_log_math generates the test statistic using the log pdf. 
  def ratio_log_math(self, sims):
    null = multivariate_normal.logpdf(sims, self.mean, self.null_cov)
    alt = multivariate_normal.logpdf(sims, self.mean, self.alt_cov)
    ratio = alt - null
    self.results = ratio

  #ratio_math generates the likelihood ratio for the importance sampling weights.
  def ratio_math(self,sims):
    null = multivariate_normal.pdf(sims, self.mean, self.null_cov)
    alt = multivariate_normal.pdf(sims, self.mean, self.alt_cov)
    ratio = alt/null
    self.results = ratio

  #set_mvalues returns the posterior marginal for each trait 
  def set_mvalues(self,sims):
    print('sims shape: ' + str(sims.shape))
    all_configs = []
    for i in range(0,len(self.alt_configs)):
      #print('i: ' + str(i))
      #print(self.alt_configs[i])
      #value = multivariate_normal.logpdf(sims, self.mean, self.alt_configs[i])
      value = multivariate_normal.pdf(sims, self.mean, self.alt_configs[i])
      #print(value[:10]) 
      if len(all_configs) is 0:
        all_configs = value
      else:
        all_configs = Likelihood_Ratio.sumlog(all_configs, value)
        #all_configs = np.add(all_configs,value)
      #for j in range(0,len(self.groups)):
        #if i in self.groups[j]:
          #if len(self.mvalues[j]) is 0:
            #self.mvalues[j] = value
          #else:
            #self.mvalues[j] = Likelihood_Ratio.sumlog(self.mvalues[j], value)
            #self.mvalues[j] = np.add(self.mvalues[j],value)
    #print(all_configs[:10])
    #self.mvalues = np.exp(np.subtract(self.mvalues,all_configs))
    #for k in range(0,len(self.mvalues)):
        #print(self.mvalues[k][:10])
        #self.mvalues[k] = np.exp(self.mvalues[k]-all_configs))
    #return self.mvalues
    return all_configs

  def get_fullprob(self,sims):
    value = multivariate_normal.pdf(sims, self.mean, self.alt_cov)
    #value = multivariate_normal.logpdf(sims, self.mean, self.alt_cov)
    return value

  #returns the results      
  def get_values(self):
    return np.asarray(self.results)

if __name__=="__main__":
    Main()
