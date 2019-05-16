'''
Author: Kodi Collins
Email: kodicollins@ucla.edu
This script takes in all input and writes output
'''

import argparse
import numpy as np


class IO():
  
  #TODO: figure out which matter
  def __init__(self):
    self.mean = []
    #model sigmaE & sigmaG
    self.sigmaE = []
    self.sigmaG = []
    #importance sampling sigmaE
    self.impSigmaE = []      
    self.combo = []
    #configurations of SigmaG             
    self.cfigSigmaG = []       
    #assign sigmaG configs to traits
    self.groups = []
    #all traits z-scores, etc.
    self.combo = []
    
    self.alpha = 1.0
    self.num = 1e6
    self.migwas = False
    self.out = ''

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
    optional.add_argument('-s', '--sampling', dest = 'sampling', default = 8, type = int,
      required = False, help = 'Variance in the the null simulations. The default is 8.')
    optional.add_argument('-x', '--num', dest = 'num', default = 1e6, type = float,
      help = 'number of simulations to run. The default is 1e6.')
    optional.add_argument('-o', '--out', dest = 'out', default = 'gwas_results.txt',
      help = 'text file where the output will be written. The default is gwas_results.txt')
    optional.add_argument('-a', '--alpha', dest = 'alpha', defalut = 1.0, type = float,
      help = 'scaling factor for mvalues')
    optional.add_argument('-m', '--migwas', action = 'store_true', default = False,
      help = 'indicates you want migwas result computed')
    args = parser.parse_args()
    read_parser(args)

  #read_parser reads the command line arguments and does some sanity checking
  def read_parser(self, args):
    env = [line.split() for line in open(args.envir)]
    gen = [line.split() for line in open(args.genetic)]
    size = [line.split() for line in open(args.pop_sizes)]
    files = [line.split() for line in open(args.gwas_files)]
    k = len(files)
    imp_var = float(args.sampling)
    self.out = args.out
    self.num = args.num
    self.alpha = args.alpha

    #few sanity checks
    if size is None or env is None or gen is None or files is None:
      raise ValueError("-envir, -genetic, -pop_sizes, and -gwas_files are required inputs, at least one is missing")
    if len(size) != len(env) and len(size) != len(gen) and len(size) != len(files):
      raise ValueError("Number of Pairwise Comparisons are not equal across input files")
    if (k*(k - 1)/2.0) != len(size):
      raise ValueError("expected number of trait pairs does not match number of lines in --pop_size file")

  #initializing attributes
  def init_self(env, gen, size, files, k, imp_var):
    #2^k models for k traits    
    include = [list(i) for i in list(itertools.product([1.0, 0.0],repeat=k))]
    #get the number of SNPs to change genetic variance to per snp              
    set_persnp()           
    
    #intermediate values (input may not be same order, etc)
    e = np.ones((k, k))
    imp = np.ones((k, k))   
    g = np.ones((k, k))   
    p = np.ones((k, k))
    #actual values
    self.mean = np.zeros(k)   
    self.sigmaG = np.ones((k, k))              
    self.sigmaE = np.ones((k, k))              
    self.impSigmaE = np.ones((k, k))               
    self.groups = [[] for i in range(k)]
    
    count = dict()          
    num = 0 
    for j in file:         
      count[j[0]] = num       
      num += 1
        
    #stores in list of lists, may not all be in same order.               
    for i in range(0,len(env)):            
      #environmental covariance           
      env_name1 = env[i][0]              
      env_name2 = env[i][1]              
      env_corr = float(env[i][2])        
      if env_corr < -1 or env_corr > 1:       
        raise ValueError("The environmenantal correlation not between -1 and 1" 
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
      gen_var1 = float(gen[i][2])
      gen_var1 = gen_var1/float(count)           
      gen_var2 = float(gen[i][3])
      gen_var2 = gen_var2/float(count)           
      gen_corr = float(gen[i][4])
      gen_cov = gen_corr*math.sqrt(gen_var1)*math.sqrt(gen_var2)
        
      one = count[gen_name1]  
      two = count[gen_name2]      
      g[one][one] = gen_var1      
      g[two][two] = gen_var2      
      g[one][two] = gen_cov       
      g[two][one] = gen_cov       
      
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
    for i in range(0,k):           
      for j in range(0,i+1):       
        if i == j:  
        self.sigmaE[i][i] = e[i][i]  
        self.sigmaG[i][i] = p[i][i]*g[i][i]
        self.impSigmaE[i][i] = imp[i][i]
        else:               
          self.sigmaE[i][j] = p[i][j]*e[i][j]
          self.sigmaE[j][i] = sigmaE[i][j]
          self.sigmaG[i][j] = math.sqrt(p[i][i])*math.sqrt(p[j][j])*g[i][j]
          self.sigmaG[j][i] = math.sqrt(p[i][i])*math.sqrt(p[j][j])*g[i][j]
        
    #setting up alt config for mvalues and assigning index to group             
    for z in range(0,len(include)):
        alt = np.copy(self.SigmaG)
        loc = include[z]     
        for i in range(0,k):   
          for j in range(0,i+1):  
            if loc[i] == 0 or loc[j] == 0:          
              alt[i,j] = 0.0
              alt[j,i] = 0.0         
        cFigSigmaG.append(alt)    
        #assign to groups   
        for l in range(0,k):           
          if loc[l] == 1.0: 
            self.groups[l].append(z)    
        
        #print('all_alt')       
        #print(alt_cov)    
        #The mean is a vector of k 0s
        mean = np.zeros(k)
        print('count: ' + str(count))
        print('divide: ' + str(divide))
        null() 

    #set_persnp gets the per SNP heritability, writes zs to file                                   
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
        gwas.to_csv(outfile, sep=' ', index=False)                                 
        del gwas


  #writes final output
  #TODO: needs to be updated to work, combo not actually created yet
  def output(self):
    if self.migwas:
      self.combo['MIGWAS_Score'] = mivalues
      self.combo['MIGWAS_Pvalue'] = mipvalues
    self.combo['LR_Score'] = lrvalues
    self.combo['LR_Pvalue'] = lrpvalues
    for i in range(len(file)):
      x = np.zeros(int(count))
      x[signif] = mvalues[i]
      self.combo['Mvalue_'+file[i][0]] = x
    self.combo.to_csv(outfile, sep=' ', index=False)
