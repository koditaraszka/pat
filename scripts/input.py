'''
Author: Kodi Collins
Email: kodicollins@ucla.edu
This script handles the parser and input
'''

import argparse
import numpy as np
import pandas
import math

class Input():
  
  #TODO: figure out which matter
  def __init__(self):
    self.mean = []
    self.thresh = 5e-8
    self.polygenic = 1e6
    #model sigmaE & sigmaG
    self.sigmaE = []
    self.sigmaG = []
    self.simSigmaG = []
    #importance sampling sigmaE
    self.impSigmaE = []      
    self.count = 1e6
    #all traits z-scores, etc.
    self.combo = []
    self.traits = []
    self.k = 0
    self.alpha = 1.0
    self.zeros = [] 
    self.migwas = False
    self.sims = False
    self.out = ''
    #self.causal = 1e6
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
    optional.add_argument('-i', '--sims', action = 'store_true', dest = 'sims', default = False, 
        help = 'for simulating data rather than analyzing --gwas_files')
    optional.add_argument('-x', '--num', dest = 'count', default = 1e6, type = float,
      help = 'number of simulations to run. The default is 1e6.')
    optional.add_argument('-y', '--polygenic', dest = 'polygenic', type = float,
      help = 'number of snps under polygenic model. If not passed in will use number of snps in real data or simulations.')
    optional.add_argument('-u', '--null', dest = 'num', default = 1e6, type = float,
      help = 'number of null simulations to run. The default is 1e6.')    
    optional.add_argument('-o', '--out', dest = 'out', default = 'gwas_results.txt',
      help = 'text file where the output will be written. The default is gwas_results.txt')
    #optional.add_argument('-a', '--alpha', dest = 'alpha', default = 1.0, type = float, help = 'TODO: depreciate')
    #optional.add_argument('-c', '--causal', dest = 'causal', default = 1e6, type = float, help = 'assuming number of causal snps, used for simulations')
    optional.add_argument('-z', '--set-zeros', dest = 'zeros', help = 'file that contains traits whose effect size is set to zero')
    optional.add_argument('-m', '--migwas', action = 'store_true', dest = 'migwas', default = False, help = 'indicates you want migwas result computed')
    optional.add_argument('-w', '--sig_thresh', dest = 'thresh', default = 5e-8, help = 'Significance threshold for M-values')
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
    #self.causal = args.causal
    imp_var = args.sampling
    self.out = args.out
    self.count = args.count
    self.num = args.num
    if args.polygenic is not None:
      self.polygenic = args.polygenic
    else:
      self.polygenic = self.count
    #self.alpha = args.alpha
    self.migwas = args.migwas
    self.thresh = float(args.thresh)
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
    sim1 = np.ones((self.k, self.k))
    sim2 = np.ones((self.k, self.k))
    sim3 = np.ones((self.k, self.k))

    #actual values
    self.mean = np.zeros(self.k)   
    self.sigmaG = np.ones((self.k, self.k))              
    self.sigmaE = np.ones((self.k, self.k))              
    self.impSigmaE = np.ones((self.k, self.k))               
    self.simSigmaG1 = np.ones((self.k, self.k))
    self.simSigmaG2 = np.ones((self.k, self.k))
    self.simSigmaG3 = np.ones((self.k, self.k))
    
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
      gen_var1 = float(gen[i][2])/float(self.polygenic)
      gen_var2 = float(gen[i][3])/float(self.polygenic)
      sim_var1_a = float(gen[i][2])/float(35000)         
      sim_var1_b = float(gen[i][2])/float(21000)         
      sim_var1_c = float(gen[i][2])/float(14000)         
      sim_var2_a = float(gen[i][3])/float(35000)     
      sim_var2_b = float(gen[i][3])/float(21000)     
      sim_var2_c = float(gen[i][3])/float(14000)     
      gen_corr = float(gen[i][4])
      gen_cov = gen_corr*math.sqrt(gen_var1)*math.sqrt(gen_var2)
      sim_cov1 = gen_corr*math.sqrt(sim_var1_a)*math.sqrt(sim_var2_a)
      sim_cov2 = gen_corr*math.sqrt(sim_var1_b)*math.sqrt(sim_var2_b)
      sim_cov3 = gen_corr*math.sqrt(sim_var1_c)*math.sqrt(sim_var2_c)
      
      one = count[gen_name1]  
      two = count[gen_name2]      
      g[one][one] = gen_var1      
      g[two][two] = gen_var2      
      g[one][two] = gen_cov       
      g[two][one] = gen_cov

      sim1[one][one] = sim_var1_a
      sim1[two][two] = sim_var2_a
      sim1[one][two] = sim_cov1
      sim1[two][one] = sim_cov1

      sim2[one][one] = sim_var1_b
      sim2[two][two] = sim_var2_b
      sim2[one][two] = sim_cov2
      sim2[two][one] = sim_cov2

      sim3[one][one] = sim_var1_c
      sim3[two][two] = sim_var2_c
      sim3[one][two] = sim_cov3
      sim3[two][one] = sim_cov3

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
          self.simSigmaG1[i][i] = p[i][i]*sim1[i][i]
          self.simSigmaG2[i][i] = p[i][i]*sim2[i][i]
          self.simSigmaG3[i][i] = p[i][i]*sim3[i][i]
        else:               
          self.sigmaE[i][j] = p[i][j]*e[i][j]
          self.sigmaE[j][i] = self.sigmaE[i][j]
          self.impSigmaE[i][j] = p[i][j]*imp[i][j]
          self.impSigmaE[j][i] = self.impSigmaE[i][j]
          self.sigmaG[i][j] = math.sqrt(p[i][i])*math.sqrt(p[j][j])*g[i][j]
          self.sigmaG[j][i] = math.sqrt(p[i][i])*math.sqrt(p[j][j])*g[i][j]
          self.simSigmaG1[i][j] = math.sqrt(p[i][i])*math.sqrt(p[j][j])*sim1[i][j]
          self.simSigmaG1[j][i] = math.sqrt(p[i][i])*math.sqrt(p[j][j])*sim1[i][j]
          self.simSigmaG2[i][j] = math.sqrt(p[i][i])*math.sqrt(p[j][j])*sim2[i][j]
          self.simSigmaG2[j][i] = math.sqrt(p[i][i])*math.sqrt(p[j][j])*sim2[i][j]
          self.simSigmaG3[i][j] = math.sqrt(p[i][i])*math.sqrt(p[j][j])*sim3[i][j]
          self.simSigmaG3[j][i] = math.sqrt(p[i][i])*math.sqrt(p[j][j])*sim3[i][j]

    if self.sims:
      print('sigmaE: ' + str(self.sigmaE))
      print('sigmaG: ' + str(self.sigmaG))
      print('simsSigmaG1: ' + str(self.simSigmaG1))
      print('simsSigmaG2: ' + str(self.simSigmaG2))
      print('simsSigmaG3: ' + str(self.simSigmaG3))
      print('sigmaE: ' + str(self.sigmaE))
      print('impSigmaE: ' + str(self.impSigmaE))
      #simCount = self.count / 22.0
      simCount1 = 5000
      simCount2 = 3000
      simCount3 = 2000
      #none
      mat = self.sigmaE
      sims = np.random.multivariate_normal(self.mean, mat, int(630000))
      print(sims.shape)
      #bdhs
      mat = self.sigmaE + self.simSigmaG1
      sims = np.append(sims,np.random.multivariate_normal(self.mean, mat, int(simCount1)), axis = 0)
      mat = self.sigmaE + self.simSigmaG2
      sims = np.append(sims,np.random.multivariate_normal(self.mean, mat, int(simCount2)), axis = 0)
      mat = self.sigmaE + self.simSigmaG3
      sims = np.append(sims,np.random.multivariate_normal(self.mean, mat, int(simCount3)), axis = 0)
      #hold1 = np.copy(self.simSigmaG1)
      hold1, hold2, hold3 = np.copy(self.simSigmaG1), np.copy(self.simSigmaG2), np.copy(self.simSigmaG3)
      #bdh
      self.simSigmaG1[3][0],self.simSigmaG1[3][1],self.simSigmaG1[3][2],self.simSigmaG1[3][3] = 0.0, 0.0, 0.0, 0.0
      self.simSigmaG1[0][3],self.simSigmaG1[1][3],self.simSigmaG1[2][3], self.simSigmaG1[3][3] = 0.0, 0.0, 0.0, 0.0
      self.simSigmaG2[3][0],self.simSigmaG2[3][1],self.simSigmaG2[3][2],self.simSigmaG2[3][3] = 0.0, 0.0, 0.0, 0.0
      self.simSigmaG2[0][3],self.simSigmaG2[1][3],self.simSigmaG2[2][3], self.simSigmaG2[3][3] = 0.0, 0.0, 0.0, 0.0
      self.simSigmaG3[3][0],self.simSigmaG3[3][1],self.simSigmaG3[3][2],self.simSigmaG3[3][3] = 0.0, 0.0, 0.0, 0.0
      self.simSigmaG3[0][3],self.simSigmaG3[1][3],self.simSigmaG3[2][3], self.simSigmaG3[3][3] = 0.0, 0.0, 0.0, 0.0
      #print('simsSigmaG1: ' + str(self.simSigmaG1))
      #print('simsSigmaG2: ' + str(self.simSigmaG2))
      #print('simsSigmaG3: ' + str(self.simSigmaG3))
      mat = self.sigmaE + self.simSigmaG1
      sims = np.append(sims,np.random.multivariate_normal(self.mean, mat, int(simCount1)), axis = 0)
      mat = self.sigmaE + self.simSigmaG2
      sims = np.append(sims,np.random.multivariate_normal(self.mean, mat, int(simCount2)), axis = 0)
      mat = self.sigmaE + self.simSigmaG3
      sims = np.append(sims,np.random.multivariate_normal(self.mean, mat, int(simCount3)), axis = 0)
      #bh
      self.simSigmaG1[1][0],self.simSigmaG1[1][1],self.simSigmaG1[1][2],self.simSigmaG1[1][3] = 0.0, 0.0, 0.0, 0.0
      self.simSigmaG1[0][1],self.simSigmaG1[1][1],self.simSigmaG1[2][1], self.simSigmaG1[3][1] = 0.0, 0.0, 0.0, 0.0
      self.simSigmaG2[1][0],self.simSigmaG2[1][1],self.simSigmaG2[1][2],self.simSigmaG2[1][3] = 0.0, 0.0, 0.0, 0.0
      self.simSigmaG2[0][1],self.simSigmaG2[1][1],self.simSigmaG2[2][1], self.simSigmaG2[3][1] = 0.0, 0.0, 0.0, 0.0
      self.simSigmaG3[1][0],self.simSigmaG3[1][1],self.simSigmaG3[1][2],self.simSigmaG3[1][3] = 0.0, 0.0, 0.0, 0.0
      self.simSigmaG3[0][1],self.simSigmaG3[1][1],self.simSigmaG3[2][1], self.simSigmaG3[3][1] = 0.0, 0.0, 0.0, 0.0
      #print('simsSigmaG1: ' + str(self.simSigmaG1))
      #print('simsSigmaG2: ' + str(self.simSigmaG2))
      #print('simsSigmaG3: ' + str(self.simSigmaG3))
      mat = self.sigmaE + self.simSigmaG1
      sims = np.append(sims,np.random.multivariate_normal(self.mean, mat, int(simCount1)), axis = 0)
      mat = self.sigmaE + self.simSigmaG2
      sims = np.append(sims,np.random.multivariate_normal(self.mean, mat, int(simCount2)), axis = 0)
      mat = self.sigmaE + self.simSigmaG3
      sims = np.append(sims,np.random.multivariate_normal(self.mean, mat, int(simCount3)), axis = 0)
      #b
      self.simSigmaG1[2][0],self.simSigmaG1[2][1],self.simSigmaG1[2][2],self.simSigmaG1[2][3] = 0.0, 0.0, 0.0, 0.0
      self.simSigmaG1[0][2],self.simSigmaG1[1][2],self.simSigmaG1[2][2], self.simSigmaG1[3][2] = 0.0, 0.0, 0.0, 0.0
      self.simSigmaG2[2][0],self.simSigmaG2[2][1],self.simSigmaG2[2][2],self.simSigmaG2[2][3] = 0.0, 0.0, 0.0, 0.0
      self.simSigmaG2[0][2],self.simSigmaG2[1][2],self.simSigmaG2[2][2], self.simSigmaG2[3][2] = 0.0, 0.0, 0.0, 0.0
      self.simSigmaG3[2][0],self.simSigmaG3[2][1],self.simSigmaG3[2][2],self.simSigmaG3[2][3] = 0.0, 0.0, 0.0, 0.0
      self.simSigmaG3[0][2],self.simSigmaG3[1][2],self.simSigmaG3[2][2], self.simSigmaG3[3][2] = 0.0, 0.0, 0.0, 0.0
      #print('simsSigmaG1: ' + str(self.simSigmaG1))
      #print('simsSigmaG2: ' + str(self.simSigmaG2))
      #print('simsSigmaG3: ' + str(self.simSigmaG3))
      mat = self.sigmaE + self.simSigmaG1
      sims = np.append(sims,np.random.multivariate_normal(self.mean, mat, int(simCount1)), axis = 0)
      mat = self.sigmaE + self.simSigmaG2
      sims = np.append(sims,np.random.multivariate_normal(self.mean, mat, int(simCount2)), axis = 0)
      mat = self.sigmaE + self.simSigmaG3
      sims = np.append(sims,np.random.multivariate_normal(self.mean, mat, int(simCount3)), axis = 0)
      #h
      self.simSigmaG1[0][0],self.simSigmaG1[0][1],self.simSigmaG1[0][2],self.simSigmaG1[0][3] = 0.0, 0.0, 0.0, 0.0
      self.simSigmaG1[0][0],self.simSigmaG1[1][0],self.simSigmaG1[2][0], self.simSigmaG1[3][0] = 0.0, 0.0, 0.0, 0.0
      self.simSigmaG2[0][0],self.simSigmaG2[0][1],self.simSigmaG2[0][2],self.simSigmaG2[0][3] = 0.0, 0.0, 0.0, 0.0
      self.simSigmaG2[0][0],self.simSigmaG2[1][0],self.simSigmaG2[2][0], self.simSigmaG2[3][0] = 0.0, 0.0, 0.0, 0.0
      self.simSigmaG3[0][0],self.simSigmaG3[0][1],self.simSigmaG3[0][2],self.simSigmaG3[0][3] = 0.0, 0.0, 0.0, 0.0
      self.simSigmaG3[0][0],self.simSigmaG3[1][0],self.simSigmaG3[2][0], self.simSigmaG3[3][0] = 0.0, 0.0, 0.0, 0.0
      self.simSigmaG1[2][2] = np.copy(hold1[2][2])
      self.simSigmaG2[2][2] = np.copy(hold2[2][2])
      self.simSigmaG3[2][2] = np.copy(hold3[2][2])
      #print('simsSigmaG1: ' + str(self.simSigmaG1))
      #print('simsSigmaG2: ' + str(self.simSigmaG2))
      #print('simsSigmaG3: ' + str(self.simSigmaG3))
      mat = self.sigmaE + self.simSigmaG1
      sims = np.append(sims,np.random.multivariate_normal(self.mean, mat, int(simCount1)), axis = 0)
      mat = self.sigmaE + self.simSigmaG2
      sims = np.append(sims,np.random.multivariate_normal(self.mean, mat, int(simCount2)), axis = 0)
      mat = self.sigmaE + self.simSigmaG3
      sims = np.append(sims,np.random.multivariate_normal(self.mean, mat, int(simCount3)), axis = 0)
      #dhs
      self.simSigmaG1 = np.copy(hold1)
      self.simSigmaG2 = np.copy(hold2)
      self.simSigmaG3 = np.copy(hold3)
      self.simSigmaG1[0][0],self.simSigmaG1[0][1],self.simSigmaG1[0][2],self.simSigmaG1[0][3] = 0.0, 0.0, 0.0, 0.0
      self.simSigmaG1[0][0],self.simSigmaG1[1][0],self.simSigmaG1[2][0], self.simSigmaG1[3][0] = 0.0, 0.0, 0.0, 0.0
      self.simSigmaG2[0][0],self.simSigmaG2[0][1],self.simSigmaG2[0][2],self.simSigmaG2[0][3] = 0.0, 0.0, 0.0, 0.0
      self.simSigmaG2[0][0],self.simSigmaG2[1][0],self.simSigmaG2[2][0], self.simSigmaG2[3][0] = 0.0, 0.0, 0.0, 0.0
      self.simSigmaG3[0][0],self.simSigmaG3[0][1],self.simSigmaG3[0][2],self.simSigmaG3[0][3] = 0.0, 0.0, 0.0, 0.0
      self.simSigmaG3[0][0],self.simSigmaG3[1][0],self.simSigmaG3[2][0], self.simSigmaG3[3][0] = 0.0, 0.0, 0.0, 0.0
      #print('simsSigmaG1: ' + str(self.simSigmaG1))
      #print('simsSigmaG2: ' + str(self.simSigmaG2))
      #print('simsSigmaG3: ' + str(self.simSigmaG3))
      mat = self.sigmaE + self.simSigmaG1
      sims = np.append(sims,np.random.multivariate_normal(self.mean, mat, int(simCount1)), axis = 0)
      mat = self.sigmaE + self.simSigmaG2
      sims = np.append(sims,np.random.multivariate_normal(self.mean, mat, int(simCount2)), axis = 0)
      mat = self.sigmaE + self.simSigmaG3
      sims = np.append(sims,np.random.multivariate_normal(self.mean, mat, int(simCount3)), axis = 0)
      #ds
      self.simSigmaG1[2][0],self.simSigmaG1[2][1],self.simSigmaG1[2][2],self.simSigmaG1[2][3] = 0.0, 0.0, 0.0, 0.0
      self.simSigmaG1[0][2],self.simSigmaG1[1][2],self.simSigmaG1[2][2], self.simSigmaG1[3][2] = 0.0, 0.0, 0.0, 0.0
      self.simSigmaG2[2][0],self.simSigmaG2[2][1],self.simSigmaG2[2][2],self.simSigmaG2[2][3] = 0.0, 0.0, 0.0, 0.0
      self.simSigmaG2[0][2],self.simSigmaG2[1][2],self.simSigmaG2[2][2], self.simSigmaG2[3][2] = 0.0, 0.0, 0.0, 0.0
      self.simSigmaG3[2][0],self.simSigmaG3[2][1],self.simSigmaG3[2][2],self.simSigmaG3[2][3] = 0.0, 0.0, 0.0, 0.0
      self.simSigmaG3[0][2],self.simSigmaG3[1][2],self.simSigmaG3[2][2], self.simSigmaG3[3][2] = 0.0, 0.0, 0.0, 0.0
      #print('simsSigmaG1: ' + str(self.simSigmaG1))
      #print('simsSigmaG2: ' + str(self.simSigmaG2))
      #print('simsSigmaG3: ' + str(self.simSigmaG3))
      mat = self.sigmaE + self.simSigmaG1
      sims = np.append(sims,np.random.multivariate_normal(self.mean, mat, int(simCount1)), axis = 0)
      mat = self.sigmaE + self.simSigmaG2
      sims = np.append(sims,np.random.multivariate_normal(self.mean, mat, int(simCount2)), axis = 0)
      mat = self.sigmaE + self.simSigmaG3
      sims = np.append(sims,np.random.multivariate_normal(self.mean, mat, int(simCount3)), axis = 0)
      print(sims.shape)
      print(self.count)
      #print(sims)
      zs = []
      bp = []
      for i in range(int(self.count)):
        bp.append([int(i)])
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
