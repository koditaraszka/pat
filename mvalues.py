from sklearn.cluster import KMeans
from scipy.stats import multivariate_normal, norm
import random
import math
import numpy as np
import itertools
import argparse
import pandas

#This class contains the pat and migwas methods
class Methods():
  def __init__(self):

  #interpretation
  def mvalues(self, zvals, loc, chrm, bp):
    k = zvals.shape[1]
    include = [list(i) for i in list(itertools.product([1.0, 0.0], repeat=k))]
    alphaData = self.prune(loc, chrm, bp)
    self.set_alpha(alphaData)
    #print(self.alpha)
    sigmaG = (self.count/self.alpha)*self.sigmaG
    mat = self.sigmaE + sigmaG
    value = multivariate_normal.logpdf(zvals, self.mean, mat)
    all_configs = np.copy(value)
    mvalues = [np.copy(value) for i in range(k)]
    for z in range(1,len(include)):
      alt = np.copy(mat)
      loc = include[z]
      for i in range(0,k):
        for j in range(0,i+1):
          if loc[i] == 0 or loc[j] == 0:
            alt[i,j] = self.sigmaE[i,j]
            alt[j,i] = self.sigmaE[j,i]
      value = multivariate_normal.logpdf(zvals, self.mean, alt)
      all_configs = Methods.sumlog(all_configs, value)
      for m in range(0,k):
        if loc[m] == 1.0:
          mvalues[m] = Methods.sumlog(mvalues[m],value)
    mvalues = np.exp(np.subtract(mvalues,all_configs))
    return mvalues

  #set alpha based on "independent" SNPs
  def prune(self, loc, chrm, bp):
    if self.sims:
      data=np.delete(loc,[chrm,bp], axis=1)
      return data
    select = []
    for i in range(1,23):
      pruningChrm = loc[np.squeeze(np.where(loc[:,chrm] == i), axis = 1)]
      maxBP = np.amax(pruningChrm, axis = 0)[bp]
      top = (maxBP+100000) - (maxBP % 100000)
      groups = [float(i) for i in range(1, int(top), 100000)]
      for j in range(0, len(groups)-1):
        lower = groups[j]
        upper = groups[j+1]-1.0
        pruningBP = pruningChrm[np.squeeze(np.where((pruningChrm[:,bp] >= lower) & (pruningChrm[:,bp] <= upper)), axis = 1)]
        if pruningBP.ndim == 1:
          select.append(np.delete(pruningBP,[chrm,bp]))
        else:
          x=pruningBP.shape[0]
          if x > 0:
            keep = random.randint(0,(x-1))
            select.append(np.delete(pruningBP[keep,:],[chrm,bp]))
      lower = groups[-1]
      upper = maxBP
      pruningBP = pruningChrm[np.squeeze(np.where((pruningChrm[:,bp] >= lower) & (pruningChrm[:,bp] <= upper)), axis = 1)]
      if pruningBP.ndim == 1:
        select.append(np.delete(pruningBP,[chrm,bp]))
      else:
        x=pruningBP.shape[0]
        if x > 0:
          keep = random.randint(0,(x-1))
          select.append(np.delete(pruningBP[keep,:],[chrm,bp]))
    select = np.array(select)
    return select

  #set alpha
  def set_alpha(self, data):
    maximum = -np.inf
    for x in range(int(self.count),0,25):
      if x == 0:
        x = 1
      x = float(x)
      mult = self.count/x
      covar = self.sigmaE + mult*self.sigmaG
      pdf = multivariate_normal.logpdf(data, self.mean, covar)
      total = np.sum(pdf)
      if total > maximum:
        maximum = total
        self.alpha = x

class Input():
  
  #TODO: figure out which matter
  def __init__(self):
    self.mean = []
    #model sigmaE & sigmaG
    self.sigmaE = []
    self.sigmaG = []
    self.simSigmaG = []
    #importance sampling sigmaE
    self.count = 1e6
    #all traits z-scores, etc.
    self.combo = []
    self.traits = []
    self.k = 0
    self.alpha = 1.0
    self.zeros = [] 
    self.out = ''
    self.causal = 1e6
    self.num = 1e6
    self.define_parser()

  #define_parser is argparser
  def define_parser(self):
    parser = argparse.ArgumentParser(description = 'This program performs multi-trait GWAS using summary statistics and interprets the results.')
  
    #required
    required = parser.add_argument_group('Required Arguments')
    required.add_argument('-e', '--envir', dest = 'envir', required = True,
      help = 'text file where each line is whitespace separated name of two traits and their pairwise environmental correlation')
    required.add_argument('-g', '--genetic', dest = 'genetic', required = True,
      help = 'text file where each line is whitespace separated name of two traits, genetic variance trait 1 and trait 2, and their genetic correlation')
    required.add_argument('-n','--pop_sizes', dest = 'pop_sizes', required = True,
      help='text file where each line is the names of two traits, sample size for trait 1 and trait 2, and the overlapping individuals')
    required.add_argument('-f', '--gwas_files', dest = 'gwas_files', required = True,
      help = 'text file where each line is the name of a trait and the path to its summary statistics')    
 
    #optional
    optional = parser.add_argument_group('Optional Arguments')
    optional.add_argument('-s', '--sampling', dest = 'sampling', default = 8, type = float,
      required = False, help = 'Variance in the the null simulations. The default is 8.')
    optional.add_argument('-x', '--num', dest = 'count', default = 1e6, type = float,
      help = 'number of simulations to run. The default is 1e6.')
    optional.add_argument('-y', '--polygenic', dest = 'polygenic', type = float,
      help = 'number of snps under polygenic model. If not passed in will use number of snps in real data or simulations.')
    optional.add_argument('-u', '--null', dest = 'num', default = 1e6, type = float,
      help = 'number of null simulations to run. The default is 1e6.')    
    optional.add_argument('-o', '--out', dest = 'out', default = 'gwas_results.txt',
      help = 'text file where the output will be written. The default is gwas_results.txt')
    #optional.add_argument('-a', '--alpha', dest = 'alpha', default = 1.0, type = float, help = 'TODO: depreciate')
    optional.add_argument('-c', '--causal', dest = 'causal', default = 1e6, type = float, help = 'assuming number of causal snps, used for simulations')
    optional.add_argument('-z', '--set-zeros', dest = 'zeros', help = 'file that contains traits whose effect size is set to zero')
    optional.add_argument('-i', '--migwas', action = 'store_true', dest = 'migwas', default = False, help = 'indicates you want migwas result computed')
    optional.add_argument('-m', '--sims', action = 'store_true', dest = 'sims', default = False, help = 'fis simulating data rather than  --gwas_files')
    args = parser.parse_args()
    self.read_parser(args)

  #read_parser reads the command line arguments and does some sanity checking
  def read_parser(self, args):
    env = [line.split() for line in open(args.envir)]
    gen = [line.split() for line in open(args.genetic)]
    size = [line.split() for line in open(args.pop_sizes)]
    files = [line.split() for line in open(args.gwas_files)]
    self.zero = []
    if args.zeros is not None:
      self.zeros = [line.split()[0] for line in open(args.zeros)]
    self.k = len(files)
    if args.sims is False and len(files) is 0:
      raise ValueError("need to pass in list of gwas files or indicate you want to simulate data")
    self.sims = args.sims
    self.causal = args.causal
    imp_var = args.sampling
    self.out = args.out
    self.count = args.count
    self.num = args.num
    #self.alpha = args.alpha
    self.migwas = args.migwas

    #few sanity checks
    if size is None or env is None or gen is None or files is None:
      raise ValueError("-envir, -genetic, -pop_sizes, and -gwas_files are required inputs, at least one is missing")
    if len(size) != len(env) and len(size) != len(gen):
      raise ValueError("Number of Pairwise Comparisons are not equal across input files")
    if (self.k*(self.k - 1)/2.0) != len(size):
      raise ValueError("expected number of trait pairs does not match number of lines in --pop_size file")
    
    #get the number of SNPs to change genetic variance to per snp              
    self.set_persnp(files)           
    self.init_self(env, gen, size, files, imp_var)
  
  #initializing attributes
  def init_self(self, env, gen, size, files, imp_var):
    
    #intermediate values (input may not be same order, etc)
    e = np.ones((self.k, self.k))
    imp = np.ones((self.k, self.k))   
    g = np.ones((self.k, self.k))   
    p = np.ones((self.k, self.k))
    sim = np.ones((self.k, self.k))

    #actual values
    self.mean = np.zeros(self.k)   
    self.sigmaG = np.ones((self.k, self.k))              
    self.sigmaE = np.ones((self.k, self.k))              
    self.impSigmaE = np.ones((self.k, self.k))               
    self.simSigmaG = np.ones((self.k, self.k))
    
    count = dict()          
    num = 0 
    for j in files:         
      count[j[0]] = num       
      num += 1
        
    #stores in list of lists, may not all be in same order.               
    for i in range(0,len(env)):            
      #environmental covariance           
      env_name1 = env[i][0]              
      env_name2 = env[i][1]              
      env_corr = float(env[i][2])        
      if env_corr < -1 or env_corr > 1:       
        raise ValueError("The environmenantal correlation not between -1 and 1") 
      one = count[env_name1]              
      two = count[env_name2]              
        
      e[one][one] = 1.0   
      e[two][two] = 1.0               
      e[one][two] = env_corr          
      e[two][one] = env_corr          
        
      #importance sampling covariance           
      imp[one][one] = imp_var  
      imp[two][two] = imp_var      
      imp[one][two] = env_corr*imp_var             
      imp[two][one] = env_corr*imp_var             
        
      #genetic covariance         
      gen_name1 = gen[i][0] 
      gen_name2 = gen[i][1]      
      gen_var1 = float(gen[i][2])/float(self.count)
      sim_var1 = float(gen[i][2])/float(self.causal)         
      gen_var2 = float(gen[i][3])/float(self.count)
      sim_var2 = float(gen[i][3])/float(self.causal)     
      gen_corr = float(gen[i][4])
      if gen_name1 in self.zeros:
        sim_var1 = 0.0
      if gen_name2 in self.zeros:
        sim_var2 = 0.0
      gen_cov = gen_corr*math.sqrt(gen_var1)*math.sqrt(gen_var2)
      sim_cov = gen_corr*math.sqrt(sim_var1)*math.sqrt(sim_var2)
      one = count[gen_name1]  
      two = count[gen_name2]      
      g[one][one] = gen_var1      
      g[two][two] = gen_var2      
      g[one][two] = gen_cov       
      g[two][one] = gen_cov       
      
      sim[one][one] = sim_var1
      sim[two][two] = sim_var2
      sim[one][two] = sim_cov
      sim[two][one] = sim_cov

      #sample size 
      pop_name1 = size[i][0] 
      pop_name2 = size[i][1]             
      pop_size1 = float(size[i][2])      
      pop_size2 = float(size[i][3])      
      pop_shared = float(size[i][4])     
      pop_ratio = pop_shared/(math.sqrt(pop_size1)*math.sqrt(pop_size2))
      one = count[pop_name1]  
      two = count[pop_name2]  
      p[one][one] = pop_size1 
      p[two][two] = pop_size2 
      p[one][two] = pop_ratio 
      p[two][one] = pop_ratio 
        
    #This sets the sigmaG and sigmaE values
    for i in range(0,self.k):           
      for j in range(0,i+1):       
        if i == j:  
          self.sigmaE[i][i] = e[i][i]  
          self.sigmaG[i][i] = p[i][i]*g[i][i]
          self.impSigmaE[i][i] = imp[i][i]
          self.simSigmaG[i][i] = p[i][i]*sim[i][i]
        else:               
          self.sigmaE[i][j] = p[i][j]*e[i][j]
          self.sigmaE[j][i] = self.sigmaE[i][j]
          self.impSigmaE[i][j] = p[i][j]*imp[i][j]
          self.impSigmaE[j][i] = self.impSigmaE[i][j]
          self.sigmaG[i][j] = math.sqrt(p[i][i])*math.sqrt(p[j][j])*g[i][j]
          self.sigmaG[j][i] = math.sqrt(p[i][i])*math.sqrt(p[j][j])*g[i][j]
          self.simSigmaG[i][j] = math.sqrt(p[i][i])*math.sqrt(p[j][j])*sim[i][j]
          self.simSigmaG[j][i] = math.sqrt(p[i][i])*math.sqrt(p[j][j])*sim[i][j]
    if self.sims:
      print('sigmaG: ' + str(self.sigmaG))
      print('simsSigmaG: ' + str(self.simSigmaG))
      print('sigmaE: ' + str(self.sigmaE))
      print('impSigmaE: ' + str(self.impSigmaE))
      #exit()
      mat = self.sigmaE + self.simSigmaG
      sims = np.random.multivariate_normal(self.mean, mat, int(self.count))
      zs = []
      bp = []
      for i in range(int(self.count)):
        bp.append([(i*100000)])
      bp = np.array(bp)
      data = np.full(int(self.count), 1)
      data = np.reshape(data, (int(self.count), -1))
      data = np.concatenate((data, bp), axis=1)
      data = np.concatenate((data, sims), axis=1)
      self.combo = pandas.DataFrame(data)
      for i in range(1,self.k+1):
        zs.append('Z_'+str(i))
      names = ['CHR','BP']
      names += zs
      self.combo.columns = names
      self.traits = [str(i) for i in range(1,self.k+1)]

  #set_persnp gets the per SNP heritability, writes zs to file
  def set_persnp(self, fileNames):
    if self.sims is False:
      traits = []
      files = []
      Z = []
      P = []
      for i in fileNames:
        traits.append(i[0])
        files.append(i[1])
        Z.append('Z_'+i[0])
        P.append('P_'+i[0])
        self.traits = traits
      gwas = pandas.read_csv(files[0], sep='\t')
      gwas = gwas.add_suffix('_'+traits[0])
      gwas = gwas.rename(index=str, columns = {'RSID_'+traits[0]: 'RSID', 'CHR_'+traits[0]: 'CHR', 'BP_'+traits[0]: 'BP', 'A1_'+traits[0]: 'A1', 'A2_'+traits[0]: 'A2'})
    
      for i in range(1,self.k):
        t = pandas.read_csv(files[i], sep='\t')
        t = t.add_suffix('_'+traits[i])
        t = t.rename(index=str, columns = {'RSID_'+traits[i]: 'RSID', 'CHR_'+traits[i]: 'CHR', 'BP_'+traits[i]: 'BP', 'A1_'+traits[i]: 'A1', 'A2_'+traits[i]: 'A2'})
        gwas = gwas.merge(t, how = 'inner', on=['RSID','CHR','BP','A1','A2'])
      self.count = float(gwas.shape[0])
      self.combo = gwas

#CHR BP Z_bmi Z_diabp Z_height Z_sysbp PAT_Pvalue
