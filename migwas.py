import numpy as np
import math

#This class computes the MI GWAS math
class Migwas():
  #def __init__(self):

  #takes in the null sims, sampling weights and percents and returns critical values
  def process_null(self, sims, weigh, percent):
    micrit = np.amax(abs(sims), axis=1) 
    miorder = micrit.argsort()
    micrit = np.stack((micrit, weigh))
    micrit = micrit[:,miorder]
    micrit[1,] = np.cumsum(micrit[1,::-1])[::-1]
    micrit[1,] = np.divide(micrit[1,],total)
    mikeep = micrit[1,].size - np.searchsorted(micrit[1,::-1], percent, side='left')
    micrit = micrit[0,]
    micrit = micrit[mikeep]
    return micrit

  #pvalues maps the real data to critical values and returns the p-value
  def pvalues(self, micrit, mivalues, percent):
    mikeep = np.searchsorted(micrit, mivalues, side='left')
    mikeep = [int(i-1) if i > 0 else int(i) for i in mikeep]
    mipvalues = percent[mikeep]
    return mipvalues

if __name__=="__main__":
  Migwas()
