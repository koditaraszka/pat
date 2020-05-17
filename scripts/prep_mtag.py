from scipy.stats import norm
import pandas as pd
import sys
import numpy as np
data = pd.read_csv(sys.argv[1], sep = " ")
#causal = sys.argv[2]
n1 = 359983 #sys.argv[3]
n2 = 340162 #sys.argv[4]
n3 = 360388 #sys.argv[5]
n4 = 340159 #sys.argv[6]

p1 = (norm.sf(abs(data["Z_1"]), 0, 1))*2
p2 = (norm.sf(abs(data["Z_2"]), 0, 1))*2
p3 = (norm.sf(abs(data["Z_3"]), 0, 1))*2
p4 = (norm.sf(abs(data["Z_4"]), 0, 1))*2

z1 = data["Z_1"]
z2 = data["Z_2"]
z3 = data["Z_3"]
z4 = data["Z_4"]

#migwas = data["MIGWAS_Pvalue"]
#mi1 = data["MIGWAS_Mvalue_1"]
#mi2 = data["MIGWAS_Mvalue_2"]
#mi3 = data["MIGWAS_Mvalue_3"]
#mi4 = data["MIGWAS_Mvalue_4"]

pat = data["PAT_Pvalue"]
pat1 = data["PAT_Mvalue_1"]
pat2 = data["PAT_Mvalue_2"]
pat3 = data["PAT_Mvalue_3"]
pat4 = data["PAT_Mvalue_4"]

chrom = data["CHR"].astype(int)
bp = (data["BP"]*10).astype(int)
A1 = np.repeat("A",data.shape[0])
A2 = np.repeat("G",data.shape[0])

n1 = np.repeat(n1, data.shape[0])
n2 = np.repeat(n2, data.shape[0])
n3 = np.repeat(n3, data.shape[0])
n4 = np.repeat(n4, data.shape[0])

MAF = np.repeat(0.05,data.shape[0])
#RSID = np.core.defchararray.add(np.repeat("RS",data.shape[0]).astype(str),bp2)

trait1 = pd.DataFrame({'MAF': MAF, 'RSID':bp,'CHR':chrom,'BP':bp,'A1':A1,'A2':A2,'N':n1,'Z':z1,'P':p1})
trait2 = pd.DataFrame({'MAF': MAF, 'RSID':bp,'CHR':chrom,'BP':bp,'A1':A1,'A2':A2,'N':n2,'Z':z2,'P':p2})
trait3 = pd.DataFrame({'MAF': MAF, 'RSID':bp,'CHR':chrom,'BP':bp,'A1':A1,'A2':A2,'N':n3,'Z':z3,'P':p3})
trait4 = pd.DataFrame({'MAF': MAF, 'RSID':bp,'CHR':chrom,'BP':bp,'A1':A1,'A2':A2,'N':n4,'Z':z4,'P':p4})

trait1.to_csv("trait1.txt", index = False, header = True, sep = "\t")
trait2.to_csv("trait2.txt", index = False, header = True, sep = "\t")
trait3.to_csv("trait3.txt", index = False, header = True, sep = "\t")
trait4.to_csv("trait4.txt", index = False, header = True, sep = "\t")
