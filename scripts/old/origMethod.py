'''
Author Kodi Collins
Email: kodicollins@ucla.edu
This software analyze GWAS summary statistics across multiple traits.
There are two methods used: Multi Trait Likelihood Ratio GWAS and Multiple Independent GWAS
'''
import numpy.linalg as linalg
import bisect

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
    self.null_cov = []
    self.count = 0
    self.size = []
    self.file = []
    self.gen = []
    self.env = []
    self.lrcrit = []
    self.percent = []
    self.micrit = []
    self.lrvalues = []
    self.mivalues = []
    self.lrpvalues = []
    self.mipvalues = []
    self.combo = []
    self.include = []
    self.nullprior = 1.0
    self.altprior = []
    self.k = 0
    self.prior_flag = None
    self.model_flag = 'full'
    self.outfile = 'gwas_results.txt'
    self.define_parser()
    print('parsed')
    self.simulate()
    print('simulated')
    self.set_percent()
    self.analyze()
    self.output()

  #define_parser is the command line argument parser. 
  #There are four required arguments -envir, -genetic, -pop_sizes, -gwas_files
  #There are three optional arguments -prior, -model, and -out
  #The help statement for each describes the argument
  def define_parser(self):
    parser = argparse.ArgumentParser(description='this program simulates a Null Panel for a Likelihood Ratio and MI GWAS to create the critical values for finding genetic variants associated with traits.  -envir, -genetic and -pop_sizes are required arguments with input text files that are whitespace separated, each line in the text file is for a set of pairwise traits. It is vital that the naming of traits in the files is consistent across files though the neither the naming choice nor the ordering have an effect.')
    parser.add_argument('-envir',help = 'text file where each line is as such trait names pairwise environmental correlation whitespace separated (eg. bmi height -0.088177 which is equivalent to height bmi -0.088177). We assume a standardizes environmental variance (ie. equal to 1). This assumption currently cannot be violated.')
    parser.add_argument('-genetic',help = 'text file where each line is the name of trait i and trait j, the genetic variance for trait i, genetic variance for trait j, and their pairwise genetic correlation (eg. bmi height 0.2468 0.453 -0.1589 which is equivalent to height bmi 0.453 0.2468 -0.1589)')
    parser.add_argument('-pop_sizes', help='text file where each line is the names of two traits and the sample size for trait 1, sample size for trait 2, and the overlapping individuals (eg. bmi height 336107 336474 336107 which is eqivalent to height bmi 336474 336107 336107')
    parser.add_argument('-gwas_files', help = 'text file where each line is the name of a trait and its corresponding GWAS summary statistics file (eg. bmi bmi.txt)')
    parser.add_argument('-prior', help = 'Will the prior be empirical (emp) or flat -prior emp or -prior flat. Flat gives all models the same weight, empirical weighs the model according to the data given. This flag does not work with -model full. If -model config is used and -prior is not used, the default is -prior flat')
    parser.add_argument('-model', help = 'Consider only the case where the  SNP affects all traits (full) or all configurations of effect (config) and weigh them according to -model. -model full or -model config, respectively. The default is -model full')
    parser.add_argument('-out', help = 'text file where the output will be written. The default is gwas_results.txt')
    parser.add_argument('-imp_samp', help = 'how much larger is the simulation')
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
    #set the other arguments if called
    if args.prior is not None:
      self.prior_flag = args.prior
    if args.model is not None:
      self.model_flag = args.model
    if args.out is not None:
      self.outfile = args.out
    if args.imp_samp is not None:
      self.imp_samp = args.imp_samp
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
    self.include = [list(i) for i in list(itertools.product([1.0, 0.0],repeat=self.k))]
    #no matter what we have to call set_prior to know how many SNPs there are to make sure the genetic variance is per SNP.
    self.set_prior()
    e = np.ones((self.k, self.k))
    g = np.ones((self.k, self.k))
    p = np.ones((self.k, self.k))
    #this matrix is the simulation matrix
    sigmaG = np.ones((self.k, self.k))
    sigmaE = np.ones((self.k, self.k))
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
      one = count[env_name1]
      two = count[env_name2]
      
      e[one][one] = 1.0
      e[two][two] = 1.0
      e[one][two] = env_corr
      e[two][one] = env_corr

      #TODO this assumes the genetic variances are for all SNPs and converts it to per SNP, how do I check this to make sure it's only doing it when it needs to
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
    #There may be a more efficiant way of setting (because this will set the pairwise twice
    #But, it seemed better to reset the values than to not set the values.
    #TODO I think if you have j go between 0 and self.k/2 or self.k/2+1 depending on if even or odd this would be sufficient
    #I have not checked this but this is not the bottleneck of our code as we only allow up to a handful of traits.
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

    #The null_cov is just the environmental effect covriance matrix
    self.null_cov = sigmaE

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
      for j in range(0,self.k):
        for k in range(0,self.k):
          alt[j,k] = sigmaG[j,k] + sigmaE[j,k]
          alt[k,j] = sigmaG[k,j] + sigmaE[k,j]
      self.alt_cov = alt
    
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
    print('null')
    #first round actually calls the classes for the first time
    #this first round sets the critical values for up 1e-1 to 1e-5 because there are enough simulations
    sims = np.random.multivariate_normal(self.mean, self.null_cov, int(1e8))
    print('round 0 done')
    #mi = Migwas(sims)
    lr = Likelihood_Ratio(sims, self.mean, self.null_cov, self.alt_cov, self.k, self.model_flag, self.nullprior, self.altprior)
    #next 99 rounds update the Migwas and Likelihood_Ratio
    for i in range(1,100):
      sims = np.random.multivariate_normal(self.mean, self.null_cov, int(1e8))
      #mi.update(sims)
      lr.update(sims)
      print('round ' + str(i) + ' done')
    #once all the data has been simulated we tell the methods and the critical values for 1e-6 to 1e-8 are set, once this is done we have all 63900 critical values returned
    #self.micrit = mi.final_set()
    self.lrcrit = lr.final_set()
    

  #set_percent begins the process of matchign the critical values to a p-value (1-alpha)
  #there are 63,900 p-values which we have just created and stored
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
    w = list(range(999,99,-1))
    w = [round(float(i)/(1e3),3) for i in w]
    x8 = [round(float(i)*1e-7,10) for i in w]    
    self.percent = x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8
    k = self.percent.index(5.0e-8)
    self.percent = self.percent[k]
    self.lrcrit = self.lrcrit[k]
    f = open(self.outfile,'w')
    f.write(' '.join(str(e) for e in [self.percent,self.lrcrit]))
    f.close()
    exit()
  #analyze begins the actual analysis of the "real" input data
  #it takes the traits and files passed in at the beginning and calls both methods on the data, this time with the real_data = True
  #currently the method assumes the data is from the UKBiobank so the column headers are according to that.
  #TODO this will need to be generalized. I may need to have a script for preprocessing the data, maybe do this when I perform LD score regression.
  def analyze(self):
      traits = []
      files = []
      for i in self.file:
        traits.append(i[0])
        files.append(i[1])
      m = Migwas([], real_data=True)
      l = Likelihood_Ratio([], self.mean, self.null_cov, self.alt_cov, self.k, self.model_flag, self.nullprior, self.altprior, real_data = True)
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
      zs = gwas[gwas.columns[gwas.columns.to_series().str.contains('Z_')]]
      #we now get the results from these Z_scores, instead of passing in simulations we pass in real data to the methods
      self.lrvalues = l.results(zs)
      self.mivalues = m.results(zs)
      #we now get their corresponding p-values
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
      elif self.prior_flag is not None:
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
      #loop over all rows of data
      for i in range(0,int(self.count)):
          #i is >= everything left of it
          x = bisect.bisect_left(self.lrcrit, self.lrvalues[i])
          x = x - 1
          #if there isn't a critical value smaller than the real data, the pvalue is 1.0
          if x == -1:
            self.lrpvalues.append(1.0)
          #else, we find the pvalue from the data
          else:
            self.lrpvalues.append(self.percent[x])
          #repeat for migwas
          y = bisect.bisect_left(self.micrit, self.mivalues[i])
          y = y - 1
          if y == -1:
            self.mipvalues.append(1.0)
          else:
            self.mipvalues.append(self.percent[y])

  #output now adds the MIGWAS and LR information the pandas table that combined all of the traits
  #it outputs it the outfile gwas_results.txt or the file name chosen. It is a whitespace separated file
  def output(self):
      self.combo['MIGWAS_Score'] = self.mivalues
      self.combo['MIGWAS_Pvalue'] = self.mipvalues
      self.combo['LR_Score'] = self.lrvalues
      self.combo['LR_Pvalue'] = self.lrpvalues
      self.combo.to_csv(self.outfile, sep=' ', index=False)



#Migwas is the class for computing the MIGWAS test statistic
#MIGWAS works by taking the largest absolute value of the set of z-scores
class Migwas:
  #MIGWAS only needs the zscores, and optionally to know if the input data is real or simulated
  def __init__(self,sims, real_data = False):
    self.sets = []
    self.crit = []
    self.hold = []
    if not real_data:
      self.critical_values(sims)
      self.lower()
      self.upper()

  #results is called on real data and gets the test statistic. It does not sort the z-scores (here called sims though not simulations) because they are real data. Once that is done, return the values
  def results(self, sims):
    self.sets = np.amax(abs(sims), axis=1)
    return self.sets

  #critical value gets the critical value and then sorts it, this is done to create the null panel from the test statistics
  def critical_values(self, sims):
    value = np.amax(abs(sims), axis=1)
    self.sets = np.sort(value)

  #update calls critical_values to get the critical values. When the class is initialized, it also calls lower to set the critical values for p-values 9.999e-1 to 1.000e-5. Upper combines the largest critical values in order to create the eventually create the critical values for p-values 9.999e-6 to 1.00e-8
  def update(self, sims):
    self.critical_values(sims)
    self.upper()

  #lower is called on the first 1e8 simulations. It collects the critical value for p-values 9.999e-1 to 1.000e-5
  def lower(self):
    one = self.sets[int(1e4):(int(9e7)+1):int(1e4)]
    one = one.tolist()
    two = self.sets[int(9.0001e7):(int(9.9e7)+1):int(1e3)]
    two = two.tolist()
    three = self.sets[int(9.90001e7):(int(9.99e7)+1):int(1e2)]
    three = three.tolist()
    four = self.sets[int(9.990001e7):(int(9.999e7)+1):int(1e1)]
    four = four.tolist()
    five = self.sets[int(9.9990001e7):(int(9.9999e7)+1):int(1e0)]
    five = five.tolist()
    #think of one = e-1, two = e-2, for which set of critical values have been generated
    self.crit = one + two + three + four + five

  #upper addresses the fact that the first round is not precise enough for the small p-values
  #we collect the largest 1e8 critical values over all simulations and store them in self.hold
  def upper(self):
    keep = np.append(self.hold, self.sets)
    keep = np.sort(keep)
    self.hold = keep[-(int(1e8)):]

  #final_set creates the final set of critical values or the critical values that correspond to 9.999e-6 to 1.00e-8
  #now that all critical values are set, return them
  def final_set(self):
    six = self.hold[int(9.990001e7):(int(9.999e7)+1):int(1e1)]
    six = six.tolist()
    seven = self.hold[int(9.9990001e7):(int(9.9999e7)+1):int(1e0)]
    seven = seven.tolist()
    eight = self.hold[int(9.9999001e7):(int(9.99999e7)+1):int(1e0)]
    eight = eight.tolist()
    #think of six = e-6, seven = e-7, etc.
    values = six + seven + eight
    #we will store the critical values in self.crit
    self.crit += values
    return self.crit


#Likelihood_Ratio is the class that computes the likelihood ratio test statistic
#Likelihood Ratio works by taking literally a likelihood ratio of the null model (only environmental effect) and the alternate model (both environmental and genetic effect). It also allows for a set of alternative models and their priors. 
class Likelihood_Ratio():
  #thie likelihood ratio method takes in a number of parameters.
  #sims is the data, where real_data says if sims are simulations or real data
  #mean, null_cov, and alt_cov are the parameters of the likelihood functions (alt_cov could be a list of alternate covariance matrices)
  #model is a flag that indicates if the model is full or config, which corresponds to if the alt_cov is a single matrix or a set of matrices, nullprior and altprior are the weights to the alt_cov if there are multiple
  def __init__(self, sims, mean, null_cov, alt_cov, k, model, nullprior, altprior, real_data = False):
    self.k = k
    self.model = model
    self.mean = mean
    self.null_cov = null_cov
    self.alt_cov = alt_cov
    self.nullprior = nullprior
    self.altprior = altprior
    self.sets = []
    self.crit = []
    self.hold = []
    if not real_data:
      self.critical_values(sims)
      self.lower()
      self.upper()

  #sumlog adds the numpy arrays of log_likelihoods. This is only needed for -model config 
  @staticmethod
  def sumlog(x, y):
    combine = np.stack((x,y), axis=1)
    high = np.amax(combine, axis=1)
    low = np.amin(combine, axis=1)
    return (np.add(low,np.log1p(np.exp(np.subtract(high,low)))))

  #results is called on the real data and gets the corresponding test statistics. It does not sort the z-scores (here called sims though they are not the null panel simulations) Once this is done it returns the values. It performs this one of two ways depending on the model_flag was set by the user (default is -model full 
  def results(self, sims):
    #this is the case where we consider all configurations of how the SNP may affect the traits
    if self.model == 'config':
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
      self.sets = ratio
    #this is the case where we are only considering the full model (the SNP affects all traits)
    else:
      null = multivariate_normal.logpdf(sims, self.mean, self.null_cov)
      alt = multivariate_normal.logpdf(sims, self.mean, self.alt_cov)
      ratio = alt - null
      self.sets = ratio
    #return the test statistics
    return self.sets

  #critical_values gets the critical values and then sorts them, this is done by creating a null panel which is passed in as sims. There are two ways of creating the critical values depending on the user specifications: consider all models or only consider the full model which is set by -model
  def critical_values(self, sims):
    #this is the case where we consider all configurations of how the SNP may affect the traits
    if self.model == 'config':
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
      self.sets = np.sort(ratio)
    #this is the case where we are only considering the full model (the SNP affects all traits)
    else:
      null = multivariate_normal.logpdf(sims, self.mean, self.null_cov)
      alt = multivariate_normal.logpdf(sims, self.mean, self.alt_cov)
      ratio = alt - null
      self.sets = np.sort(ratio)

  #update calls critical_values to get the critical values. When the class is initialized, it also calls lower to set the critical values for p-values 9.999e-1 to 1.000e-5. Upper combines the largest critical values in order to create the eventually create the critical values for p-values 9.999e-6 to 1.00e-8
  def update(self, sims):
    self.critical_values(sims)
    self.upper()

  #lower is called on the first 1e8 simulations. It collects the critical value for p-values 9.999e-1 to 1.000e-5  
  def lower(self):
    one = self.sets[int(1e4):(int(9e7)+1):int(1e4)]
    one = one.tolist()
    two = self.sets[int(9.0001e7):(int(9.9e7)+1):int(1e3)]
    two = two.tolist()
    three = self.sets[int(9.90001e7):(int(9.99e7)+1):int(1e2)]
    three = three.tolist()
    four = self.sets[int(9.990001e7):(int(9.999e7)+1):int(1e1)]
    four = four.tolist()
    five = self.sets[int(9.9990001e7):(int(9.9999e7)+1):int(1e0)]
    five = five.tolist()
    self.crit = one + two + three + four + five

  #upper addresses the fact that the first round is not precise enough for the small p-values
  #we collect the largest 1e8 critical values over all simulations and store them in self.hold
  def upper(self):
    keep = np.append(self.hold, self.sets)
    keep = np.sort(keep)
    self.hold = keep[-(int(1e8)):]

  #final_set creates the final set of critical values or the critical values that correspond to 9.999e-6 to 1.00e-8
  #now that all critical values are set, return them
  def final_set(self):
    six = self.hold[int(9.990001e7):(int(9.999e7)+1):int(1e1)]
    six = six.tolist()
    seven = self.hold[int(9.9990001e7):(int(9.9999e7)+1):int(1e0)]
    seven = seven.tolist()
    eight = self.hold[int(9.9999001e7):(int(9.99999e7)+1):int(1e0)]
    eight = eight.tolist()
    values = six + seven + eight
    self.crit += values
    return self.crit

if __name__=="__main__":
    Main()
