'''
Author Kodi Collins
Email: kodicollins@ucla.edu
This script performs GWAS by leveraging genetic correlation between traits.
Run python theMethod.py -h or see the readme for more information on how to use the method.  
There is some hardcoding in this method:
  2. the number of simulations performed is set to 1e6, searching for 'sim =' will get the right line if you want to change it
  3. the trait/environmental variance-covariance is assumed to be already standardized 
'''
import numpy.linalg as linalg
import pandas
from scipy.stats import multivariate_normal
import numpy as np
import argparse
import math
import itertools

#This is the main class. This class is responsible for input/output and for calling the two classes Migwas and Likelihood_Ratio.
#Those methods are responsible for the actual computation work.
class Main():
  def __init__(self):
    self.mean = []
    self.alt_cov = []
    self.alt_cov2 = []
    self.null_cov = []
    self.null_cov2 = []
    self.sim_cov = []
    self.sim_cov2 = []
    self.imp_var = 1.0
    self.percent = []
    self.count = 0
    self.size = []
    self.file = []
    self.gen = []
    self.env = []
    self.lrcrit = []
    self.lrcrit2 = []
    self.lrvalues = []
    self.lrvalues2 = []
    self.lrpvalues = []
    self.lrpvalues2 = []
    self.micrit = []
    self.micrit2 = []
    self.mivalues = []
    self.mivalues2 = []
    self.mipvalues = []
    self.mipvalues2 = []
    self.combo = []
    self.include = []
    self.nullprior = 1.0
    self.altprior = []
    self.k = 0
    self.prior_flag = None
    self.micount = []
    self.micount2 = []
    self.lrcount = []
    self.lrcount2 = []
    self.multcount = []
    self.model_flag = 'full'
    self.outfile = 'gwas_results.txt'
    self.zero = []
    self.alt_sim = []
    self.define_parser()
    for i in range(0,201,1):
      self.times = float(i)/200.0
      print('times')
      print(self.times)
      self.multcount.append(self.times)
      self.simulate()
      self.set_percent()
      self.analyze()
      self.output()

  #define_parser is the command line argument parser. 
  #There are four required arguments -envir, -genetic, -pop_sizes, -gwas_files
  #There are three optional arguments -prior, -model, and -out
  #The help statement for each describes the argument
  def define_parser(self):
    #print('define_parser')
    parser = argparse.ArgumentParser(description='this program simulates a Null Panel for a Likelihood Ratio and MI GWAS to create the critical values for finding genetic variants associated with traits.  -envir, -genetic and -pop_sizes are required arguments with input text files that are whitespace separated, each line in the text file is for a set of pairwise traits. It is vital that the naming of traits in the files is consistent across files though the neither the naming choice nor the ordering have an effect.')
    parser.add_argument('-envir',help = 'text file where each line is as such trait names pairwise environmental correlation whitespace separated (eg. bmi height -0.088177 which is equivalent to height bmi -0.088177). We assume a standardizes environmental variance (ie. equal to 1). This assumption currently cannot be violated.')
    parser.add_argument('-genetic',help = 'text file where each line is the name of trait i and trait j, the genetic variance for trait i, genetic variance for trait j, and their pairwise genetic correlation (eg. bmi height 0.2468 0.453 -0.1589 which is equivalent to height bmi 0.453 0.2468 -0.1589)')
    parser.add_argument('-pop_sizes', help='text file where each line is the names of two traits and the sample size for trait 1, sample size for trait 2, and the overlapping individuals (eg. bmi height 336107 336474 336107 which is eqivalent to height bmi 336474 336107 336107')
    parser.add_argument('-gwas_files', help = 'text file where each line is the name of a trait and its corresponding GWAS summary statistics file (eg. bmi bmi.txt)')
    parser.add_argument('-prior', help = 'Will the prior be empirical (emp) or flat -prior emp or -prior flat. Flat gives all models the same weight, empirical weighs the model according to the data given. This flag does not work with -model full. If -model config is used and -prior is not used, the default is -prior flat')
    parser.add_argument('-model', help = 'Consider only the case where the  SNP affects all traits (full) or all configurations of effect (config) and weigh them according to -model. -model full or -model config, respectively. The default is -model full')
    parser.add_argument('-sampling', help = 'This takes in an number, which will be used as the variance in the the null simulations')
    parser.add_argument('-out', help = 'text file where the output will be written. The default is gwas_results.txt')
    parser.add_argument('-zero', help = 'list the traits whose variance will be set to zero. (optional)')
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
    self.imp_var = float(args.sampling)
    #set the other arguments if called
    if args.prior is not None:
      self.prior_flag = args.prior
    if args.model is not None:
      self.model_flag = args.model
    if args.out is not None:
      self.outfile = args.out
    if args.zero is not None:
      self.zero = [line.split() for line in open(args.zero)]
      self.zero = self.zero[0]
      print('which are zero')
      print(self.zero)
    #check the arguments now that they have been set
    if self.size is None or self.env is None or self.gen is None or self.file is None:
      raise ValueError("-envir, -genetic, -pop_sizes, and -gwas_files are required inputs, at least one is missing")
    if len(self.size) != len(self.env) and len(self.size) != len(self.gen) and len(self.size) != len(self.file):
      raise ValueError("Number of Pairwise Comparisons are not equal across input files")
    if (self.k*(self.k - 1)/2.0) != len(self.size):
      raise ValueError("k(k-1)/2 for the k traits from gwas_files does not equal the number of lines of pairwise comparisons in the -pop_size input file")
    if self.model_flag == 'full' and self.prior_flag is not None:
      raise ValueError("-model full is incompatible with a prior. Do not set -prior with -model full. You may set a prior with -model config")
    if self.model_flag == 'config' and self.prior_flag is None:
      self.prior_flag = 'flat'

  #simulate prepares the input data for simulation
  #it actually creates the environmental and genetic covariance matrices needed
  #the trait names should be consistent across input files (check this?)
  #simulate uses the trait names to store information in dictionaries and lists of lists
  def simulate(self):
    #if the model is config there are 2^k - 1 alternative models to consider (k is the number of traits)
    #self.include finds all 2^k, and we use if 1.0 in i to skip the model with all 0.0s which would equal self.null_cov
    if self.prior_flag is not None:
      self.include = [list(i) for i in list(itertools.product([1.0, 0.0],repeat=self.k))]
    #no matter what we have to call set_prior to know how many SNPs there are to make sure the genetic variance is per SNP.
    self.set_prior()
    e = np.ones((self.k, self.k))
    e2 = np.ones((self.k, self.k))
    g = np.ones((self.k, self.k))
    p = np.ones((self.k, self.k))
    sigmaG = np.ones((self.k, self.k))
    sigmaE = np.ones((self.k, self.k))
    sigmaE2 = np.ones((self.k, self.k))
    sim_v = np.ones((self.k, self.k))
    sigmaS = np.ones((self.k, self.k))
    sigmaS2 = np.ones((self.k, self.k))
    sim_v2 = np.ones((self.k, self.k))
    count = dict()
    num = 0
    
    for j in self.file:
      count[j[0]] = num
      num += 1
    
    #this puts the data into the list of lists of values.
    #e[trait1][trait2] == e[trait2][trait1]
    #this is true for all lists of lists
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
      e[one][two] = self.times*env_corr
      e[two][one] = self.times*env_corr
    
      e2[one][one] = 1.0
      e2[two][two] = 1.0
      e2[one][two] = 0.0
      e2[two][one] = 0.0
      #simulated covariance
      sim_v[one][one] = self.imp_var
      sim_v[two][two] = self.imp_var
      sim_v[one][two] = (self.times*env_corr)*self.imp_var
      sim_v[two][one] = (self.times*env_corr)*self.imp_var

      sim_v2[one][one] = self.imp_var
      sim_v2[two][two] = self.imp_var
      sim_v2[one][two] = 0.0*self.imp_var
      sim_v2[two][one] = 0.0*self.imp_var
      #genetic covariance
      gen_name1 = self.gen[i][0]
      gen_name2 = self.gen[i][1]
      #TODO how to actually do this and check it
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
    #There may be a more efficiant way of setting (because this will set the pairwise twice
    #But, it seemed better to reset the values than to not set the values.
    #TODO I think if you have j go between 0 and self.k/2 or self.k/2+1 depending on if even or odd this would be sufficient
    #I have not checked this but this is not the bottleneck of our code as we only allow up to a handful of traits.
    for i in range(0,self.k):
        for j in range(0,self.k):
          if i == j:
            sigmaE[i][i] = e[i][i]
            sigmaS[i][i] = sim_v[i][i]
            sigmaE2[i][i] = e2[i][i]
            sigmaS2[i][i] = sim_v2[i][i]
            sigmaG[i][i] = p[i][i]*g[i][i]
          else:
            sigmaE[i][j] = p[i][j]*e[i][j]
            sigmaS[i][j] = p[i][j]*sim_v[i][j]
            sigmaE2[i][j] = p[i][j]*e2[i][j]
            sigmaS2[i][j] = p[i][j]*sim_v2[i][j]
            sigmaE[j][i] = sigmaE[i][j]
            sigmaS[j][i] = sigmaS[i][j]
            sigmaE2[j][i] = sigmaE2[i][j]
            sigmaS2[j][i] = sigmaS2[i][j]
            sigmaG[i][j] = math.sqrt(p[i][i])*math.sqrt(p[j][j])*g[i][j]
            sigmaG[j][i] = sigmaG[i][j]

    #null_cov is just the environmental effect covariance matrix
    self.null_cov = sigmaE
    self.null_cov2 = sigmaE2
    self.alt_sim = sigmaE
    #sim_cov is the actual covariance matrix used in the simulations, if -sampling isn't used then it equals null_cov
    self.sim_cov = sigmaS
    self.sim_cov2 = sigmaS2
    #if there is a prior flag then the model is config which means there are 2^k - 1 alternative models to consider (k is the number of traits)
    #self.include finds all 2^k, and we use if 1.0 in i to skip the model with all 0.0s which would equal self.null_cov
    if self.prior_flag is not None:
      for z in range(0,len(self.include)):
        alt = np.zeros((self.k,self.k))
        i = self.include[z]
        if 1.0 in i:
          for j in range(0,self.k):
            for k in range(0,self.k):
              alt[j,k] = (i[j]*sigmaG[j,k]) + sigmaE[j,k]
              alt[k,j] = (i[j]*sigmaG[k,j]) + sigmaE[k,j]
          self.alt_cov.append(alt)
    #if there isn't a prior, we are only considering the full model which means only the alternative model with all 1s.
    else:
      if self.model_flag != 'full':
        raise ValueError('The prior flag is None and the model flag is not full. This should have been caught but was not.') 
      alt = np.zeros((self.k,self.k))
      alt2 = np.zeros((self.k,self.k))
      for j in range(0,self.k):
        for k in range(0,self.k):
          alt[j,k] = 500.0*sigmaG[j,k] + sigmaE[j,k]
          alt[k,j] = 500.0*sigmaG[k,j] + sigmaE[k,j]
          alt2[j,k] = 500.0*sigmaG[j,k] + sigmaE2[j,k]
          alt2[k,j] = 500.0*sigmaG[k,j] + sigmaE2[k,j]
      self.alt_cov = alt
      #self.alt_sim = alt
      self.alt_cov2 = alt2
      print('sim_cov')
      print(self.sim_cov)
      print('sim_cov2')
      print(self.sim_cov2)
      print('sigmaG')
      print(sigmaG)
      print('sigmaE')
      print(sigmaE)
      print('sigmaE2')
      print(sigmaE2)
      print('self.alt_cov')
      print(self.alt_cov)
      print('self.alt_cov2')
      print(self.alt_cov2)
      print('self.alt_sim')
      print(self.alt_sim)
    #The mean is a vector of k 0s
    self.mean = np.zeros(self.k)
    #we will now actually simulate the null panel
    self.null()


  #null actually simulates the data
  #we simulate here and pass into the classes Migwas and Likelihood_Ratio so that both methods have the same simulated data
  #this reduces the cost of simulating multiple times and prevents random differences in simulations from affecting our method comparison
  #we simulate in 1e8 levels due to limitations in scale and repeat for 100 iterations for a total of 1e10 simulations.
  #This allows us to report p-values from 9.999e-2 to 1.00e-8 all all but 1e-8 we have p-values for 4 significant digits and for 1e-8 "level" we have 3 digits
  #The cost difference between 100 and 1000 simulations was too high to have 4 digits for 9.99e-8 to 1.00e-8 p-values (meaning 9.999e-8 to 1.000e-8)
  def null(self):
    sims = np.random.multivariate_normal(self.mean, self.sim_cov, int(1e6))
    w = Likelihood_Ratio(sims, self.mean, self.sim_cov, self.null_cov, 'full', [], [], log=False) 
    weigh = w.get_values()
    total = np.sum(weigh)
    lr = Likelihood_Ratio(sims, self.mean, self.null_cov, self.alt_cov, self.model_flag, self.nullprior, self.altprior)
    mi = Migwas(sims)
    del sims
    self.lrcrit = lr.get_values()
    self.micrit = mi.get_values()
    lrorder = self.lrcrit.argsort()
    miorder = self.micrit.argsort()
    self.lrcrit = np.stack((self.lrcrit, weigh))
    self.micrit = np.stack((self.micrit, weigh))
    del weigh
    self.lrcrit = self.lrcrit[:,lrorder]
    self.micrit = self.micrit[:,miorder]
    self.lrcrit[1,] = np.cumsum(self.lrcrit[1,::-1])[::-1]
    self.micrit[1,] = np.cumsum(self.micrit[1,::-1])[::-1]
    del lrorder
    del miorder
    self.lrcrit[1,] = np.divide(self.lrcrit[1,],total)
    self.micrit[1,] = np.divide(self.micrit[1,],total)
    
    sims2 = np.random.multivariate_normal(self.mean, self.sim_cov2, int(1e6))
    w2 = Likelihood_Ratio(sims2, self.mean, self.sim_cov2, self.null_cov2, 'full', [], [], log=False)
    weigh2 = w2.get_values()
    total2 = np.sum(weigh2) 
    lr2 = Likelihood_Ratio(sims2, self.mean, self.null_cov2, self.alt_cov2, self.model_flag, self.nullprior, self.altprior)
    mi2 = Migwas(sims2)
    del sims2
    self.lrcrit2 = lr2.get_values()
    self.micrit2 = mi2.get_values()
    lrorder2 = self.lrcrit2.argsort()
    miorder2 = self.micrit2.argsort()
    self.lrcrit2 = np.stack((self.lrcrit2, weigh2))
    self.micrit2 = np.stack((self.micrit2, weigh2))
    del weigh2
    self.lrcrit2 = self.lrcrit2[:,lrorder2]
    self.micrit2 = self.micrit2[:,miorder2]
    self.lrcrit2[1,] = np.cumsum(self.lrcrit2[1,::-1])[::-1]
    self.micrit2[1,] = np.cumsum(self.micrit2[1,::-1])[::-1]
    del lrorder2
    del miorder2
    self.lrcrit2[1,] = np.divide(self.lrcrit2[1,],total2)
    self.micrit2[1,] = np.divide(self.micrit2[1,],total2)

  #set_percent begins the process of matchign the critical values to a p-value (1-alpha)
  #there are 72k p-values which we have just created and stored
  def set_percent(self):
    #this gets the critical vaules for each level
    #this is so we have 100 at each level as decimals and in the right order
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
    #this adds 1 to the beginning so any values smaller than the smallest critical value will map to 1.0
    #this also means that any value larger than our largest critical value will map to 1e-8
    percent = [1.0] + x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8
    self.percent = np.asarray(percent)
    self.map_values()
  
  #map values take the sum(weigh[i:])/total that is >= to the percent of interest
  def map_values(self):
    #self.percent is stored largest to smallest
    #this reverses the lrcrit[1,::-1]
    #this reverses the index self.lrcrit.size - np.searchsorted()
    #while side='right' is self.lrcrit[1,i-1] <= self.percent < self.lrcrit[1,i]
    #as we're doing it in reverse order, we want the side='left'
    #the critical value will therefore map to a pvalue that is strictly greater than its pvalue 
    lrkeep = self.lrcrit[1,].size - np.searchsorted(self.lrcrit[1,::-1],self.percent,side='left')
    lrkeep2 = self.lrcrit2[1,].size - np.searchsorted(self.lrcrit2[1,::-1],self.percent,side='left')
    mikeep = self.micrit[1,].size - np.searchsorted(self.micrit[1,::-1],self.percent,side='left')
    mikeep2 = self.micrit2[1,].size - np.searchsorted(self.micrit2[1,::-1],self.percent,side='left')
    #only keep the critical vaules
    self.lrcrit = self.lrcrit[0,]
    self.lrcrit2 = self.lrcrit2[0,]
    self.micrit = self.micrit[0,]
    self.micrit2 = self.micrit2[0,]
    #only keep the critical values that mapped
    self.lrcrit = self.lrcrit[lrkeep]
    self.lrcrit2 = self.lrcrit2[lrkeep2]
    self.micrit = self.micrit[mikeep]
    self.micrit2 = self.micrit2[mikeep2]
 
  
  #analyze begins the actual analysis of the "real" input data
  #it takes the traits and files passed in at the beginning and calls both methods on the data, this time with the real_data = True
  #currently the method assumes the data is from the UKBiobank so the column headers are according to that.
  #TODO this will need to be generalized. I may need to have a script for preprocessing the data, maybe do this when I perform LD score regression.
  def analyze(self):
    '''traits = []
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
    #self.combo is the single pandas table that has all the traits together
    self.combo = gwas
    #this is the Z_scores or s statistic from LD score regression on the genetic effects
    zs = gwas[gwas.columns[gwas.columns.to_series().str.contains('Z_')]]'''
    #we now get the results from these Z_scores, instead of passing in simulations we pass in real data to the methods
    zs = np.random.multivariate_normal(self.mean, self.alt_sim, int(1e8))
    #self.combo = pandas.DataFrame({'Z_1': zs[:,0]})
    #self.combo['Z_2'] = zs[:,1]
    #self.combo['Z_3'] = zs[:,2]
    #self.combo['Z_4'] = zs[:,3]
    m = Migwas(zs)
    self.mivalues = m.get_values()
    l = Likelihood_Ratio(zs, self.mean, self.null_cov, self.alt_cov, self.model_flag, self.nullprior, self.altprior)
    l2 = Likelihood_Ratio(zs, self.mean, self.null_cov2, self.alt_cov2, self.model_flag, self.nullprior, self.altprior)
    self.lrvalues = l.get_values()
    self.lrvalues2 = l2.get_values()
    #we now get their corresponding p-values
    del zs
    self.pvalues()

  #set_prior sets the prior if -model is config, it's still called if the model is full because we need to know the number of SNPs.
  #there are two possible priors flat and empirical
  #flat prior is going to equally way all models and empirical prior is going to weigh them according to the real data
  def set_prior(self):
    #first find out how many SNPs there are, also this is used by the empirical prior
    self.count = 7951094
    '''traits = []
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
    zs = gwas[gwas.columns[gwas.columns.to_series().str.contains('Z_|P_')]]
    del gwas

    #a flat prior is really straight forward
    if self.prior_flag == 'flat':
      flat = math.log(1.0/float(len(self.include)))
      self.nullprior = flat
      self.altprior = [flat]*(len(self.include)-1)
  
    #this is done instead of checking if == emp or empirical to forgo making sure the user calls it exactly emp
    elif self.prior_flag is 'emp':
      #we find the p-values <= 1e-4 matches this specific configuration where the remaining p-values need to be > 1e-4
      #so if there are 3 traits, bmi, diastolic blood pressure, and height and I want to see the model there only being a genetic effect on bmi and height
      #I want to find the p-values from the single trait gwas where it's <= 1e-4 for bmi and height and > 1e-4 for diastolic blood pressure
      #The null model's prior will be the case where all p-values from the single trait gwas are > 1e-4, it's the highest prior.
      reset = []
      for i in range(0,len(self.include)):
        sub = zs
        if 1.0 in self.include[i]:
          for j in range(0,self.k):
            if self.include[i][j] == 1.0:
              sub = sub.loc[sub[P[j]] <= 1e-4]
            else:
              sub = sub.loc[sub[P[j]] > 1e-4]
          x = float(sub.shape[0])/self.count
          if self.include[i].count(1.0) == 1:
            reset.append(i)
          self.altprior.append(x)
          self.nullprior -= x
        
      #average all univariate effects
      summ = 0.0
      for x in reset:
        summ += self.altprior[x]
      summ = summ/float(len(reset))
        
      #set all univariate priors equal in case only one trait has a lot of low p-values.
      for x in reset:
        self.altprior[x] = summ
        
      if self.nullprior <= 0.0:
        raise ValueError('the null prior is <= 0.0, this makes NO SENSE!')
      self.nullprior = math.log(self.nullprior)
      for i in range(0,len(self.altprior)):
        if self.altprior[i] == 0.0:
          #smallest log without error
          self.altprior[i] = math.log(1e-323)
        else:
          self.altprior[i] = math.log(self.altprior[i])'''
  
  #pvalues maps the real data to the corresponding critical value and it's pvalue
  def pvalues(self):
    #using side='left so self.micrit[i-1] <= self.mivalues < self.micrit[i]
    mikeep = np.searchsorted(self.micrit,self.mivalues, side='left')
    mikeep2 = np.searchsorted(self.micrit2,self.mivalues, side='left')
    lrkeep = np.searchsorted(self.lrcrit,self.lrvalues, side='left')
    lrkeep2 = np.searchsorted(self.lrcrit2,self.lrvalues2, side='left')
    mikeep[mikeep >= len(self.percent)] = len(self.percent) - 1
    mikeep2[mikeep2 >= len(self.percent)] = len(self.percent) - 1
    lrkeep[lrkeep >= len(self.percent)] = len(self.percent) - 1
    lrkeep2[lrkeep2 >= len(self.percent)] = len(self.percent) - 1
    self.mipvalues = self.percent[mikeep]
    self.mipvalues2 = self.percent[mikeep2]
    self.lrpvalues = self.percent[lrkeep]
    self.lrpvalues2 = self.percent[lrkeep2]
    print('count of migwas pvalues <= 5e-8')
    x1 = len(self.mipvalues[self.mipvalues <= 5e-8])
    print(x1)
    self.micount.append(x1)
    print('count of lr pvalues <= 5e-8')
    x2 = len(self.lrpvalues[self.lrpvalues <= 5e-8])
    print(x2)
    self.lrcount.append(x2)

    print('count of bad migwas pvalues <= 5e-8')
    x1 = len(self.mipvalues2[self.mipvalues2 <= 5e-8])
    print(x1)
    self.micount2.append(x1)
    print('count of bad lr pvalues <= 5e-8')
    x2 = len(self.lrpvalues2[self.lrpvalues2 <= 5e-8])
    print(x2)
    self.lrcount2.append(x2)
  #output now adds the MIGWAS and LR information the pandas table that combined all of the traits
  #it outputs it the outfile gwas_results.txt or the file name chosen. It is a whitespace separated file
  def output(self):
    self.combo = pandas.DataFrame({'multiplier': self.multcount})  
    self.combo['GoodMIGWAS_FP'] = self.micount
    self.combo['BadMIGWAS_FP'] = self.micount2
    self.combo['GoodLR_FP'] = self.lrcount
    self.combo['BadLR_FP'] = self.lrcount2
    self.combo.to_csv(self.outfile, sep=' ', index=False)

#Migwas is the class for computing the MIGWAS test statistic
#MIGWAS works by taking the largest absolute value of the set of z-scores
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
#Likelihood Ratio works by taking literally a likelihood ratio of the null model (only environmental effect) and the alternate model (both environmental and genetic effect). It also allows for a set of alternative models and their priors. 
class Likelihood_Ratio():
  #thie likelihood ratio method takes in a number of parameters.
  #sims is the data, where real_data says if sims are simulations or real data
  #mean, null_cov, and alt_cov are the parameters of the likelihood functions (alt_cov could be a list of alternate covariance matrices)
  #model is a flag that indicates if the model is full or config, which corresponds to if the alt_cov is a single matrix or a set of matrices, nullprior and altprior are the weights to the alt_cov if there are multiple
  def __init__(self, sims, mean, null_cov, alt_cov, model, nullprior, altprior, log=True):
    self.model_flag = model
    self.mean = mean
    self.null_cov = null_cov
    self.alt_cov = alt_cov
    self.nullprior = nullprior
    self.altprior = altprior
    self.results = []
    if log is False:
      self.ratio_math(sims)
    else:
      self.ratio_log_math(sims)

  #sumlog adds the numpy arrays of log_likelihoods. This is only needed for -model config 
  @staticmethod
  def sumlog(x, y):
    combine = np.stack((x,y), axis=1)
    high = np.amax(combine, axis=1)
    low = np.amin(combine, axis=1)
    return (np.add(low,np.log1p(np.exp(np.subtract(high,low)))))

  #ratio_log_math generates the test statistic using the log pdf. This is used except for when generating the weights 
  def ratio_log_math(self, sims):
    #this is the case where we consider all configurations of how the SNP may affect the traits
    if self.model_flag == 'config':
      #priors are already in terms of log
      null = multivariate_normal.logpdf(sims, self.mean, self.null_cov)
      null = null - self.nullprior
      one = multivariate_normal.logpdf(sims, self.mean, self.alt_cov[0])
      one += self.altprior[0]
      two = multivariate_normal.logpdf(sims, self.mean, self.alt_cov[1])
      two += self.altprior[1]
      alt = Likelihood_Ratio.sumlog(one, two)
      for i in range(2,len(self.alt_cov)):
        add = multivariate_normal.logpdf(sims, self.mean, self.alt_cov[i])
        add += self.altprior[i]
        alt = Likelihood_Ratio.sumlog(add, alt)
      alt = math.log(len(self.alt_cov)) + alt
      ratio = alt - null
      self.results = ratio
    #this is the case where we are only considering the full model (the SNP affects all traits)
    else:
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

  #get_values returns the test statistics, they are not sorted
  def get_values(self):
    return np.asarray(self.results)

if __name__=="__main__":
    Main()
