import pandas as pd
import sys

#k= sys.argv[1]
#causal = sys.argv[2]
direc=sys.argv[1]
pat_dir=sys.argv[2]
hipo_dir=sys.argv[3]
pat = pd.read_csv(pat_dir, sep = ' ')
hipo = pd.read_csv(hipo_dir, sep = '\t')
m1 = pd.read_csv(direc+"/mtag_results_trait_1.txt", sep = '\t')
m2 = pd.read_csv(direc+"/mtag_results_trait_2.txt", sep = '\t')
m3 = pd.read_csv(direc+"/mtag_results_trait_3.txt", sep = '\t')
m4 = pd.read_csv(direc+"/mtag_results_trait_4.txt", sep = '\t')
migwas = pat

if(pat.shape[0] != m1.shape[0]):
	print("shape does not match")

size = int((pat.shape[0] - 630000))

m1.rename(index=str, columns = {"Z": "Z1", "N":"N1","mtag_beta":"mtag_beta1", "mtag_se":"mtag_se1", "mtag_z":"mtag_z1", "mtag_pval":"mtag_pval1"}, inplace = True)
m2.rename(index=str, columns = {"Z": "Z2", "N":"N2","mtag_beta":"mtag_beta2", "mtag_se":"mtag_se2", "mtag_z":"mtag_z2", "mtag_pval":"mtag_pval2"}, inplace = True)
m3.rename(index=str, columns = {"Z": "Z3", "N":"N3","mtag_beta":"mtag_beta3", "mtag_se":"mtag_se3", "mtag_z":"mtag_z3", "mtag_pval":"mtag_pval3"}, inplace = True)
m4.rename(index=str, columns = {"Z": "Z4", "N":"N4","mtag_beta":"mtag_beta4", "mtag_se":"mtag_se4", "mtag_z":"mtag_z4", "mtag_pval":"mtag_pval4"}, inplace = True)
mtag = m1.merge(m2, how = 'inner', on = ['SNP','CHR','BP','A1','A2','FRQ'])
mtag = mtag.merge(m3, how = 'inner', on = ['SNP','CHR','BP','A1','A2','FRQ'])
mtag = mtag.merge(m4, how = 'inner', on = ['SNP','CHR','BP','A1','A2','FRQ'])
size_list = [630000,
             5000, 3000, 2000,
             5000, 3000, 2000,
             5000, 3000, 2000,
             5000, 3000, 2000,
             5000, 3000, 2000,
             5000, 3000, 2000,
             5000, 3000, 2000]
base = 0
top = base + size_list[0]-1
pat_list = [pat.loc[base:top,:]]
migwas_list = [migwas.loc[base:top,:]]
mtag_list = [mtag.loc[base:top,:]]
hipo_list = [hipo.loc[base:top,:]]
for i in range(1,22):
  base = top + 1
  top = base + size_list[i]-1
  pat_list += [pat.loc[base:top,:]]
  migwas_list += [migwas.loc[base:top,:]]
  mtag_list += [mtag.loc[base:top,:]]
  hipo_list += [hipo.loc[base:top,:]]

effect = ['none',
'bdhs','bdhs','bdhs',
'bdh','bdh','bdh',
'bh','bh','bh',
'b','b','b',
'h','h','h',
'dhs','dhs','dhs',
'ds','ds','ds']
causal_snps = ['0',
'35,000','21,000','14,000',
'35,000','21,000','14,000',
'35,000','21,000','14,000',
'35,000','21,000','14,000',
'35,000','21,000','14,000',
'35,000','21,000','14,000',
'35,000','21,000','14,000']
#50,000 25,000, 5,000
#print(len(pat))
#print(len(pat_list))
#print(pat_list)
#exit()
output = ['effect', 'causal', 'pat_sig', 'pat_b1', 'pat_d1', 'pat_h1', 'pat_s1', 'pat_b2', 'pat_d2', 'pat_h2', 'pat_s2', 'migwas_sig', 'migwas_b1', 'migwas_d1', 'migwas_h1', 'migwas_s1', 'migwas_b2', 'migwas_d2', 'migwas_h2', 'migwas_s2', 'hipo_sig', 'hipo_b1', 'hipo_d1', 'hipo_h1', 'hipo_s1', 'hipo_b2', 'hipo_d2', 'hipo_h2', 'hipo_s2', 'mtag_sig', 'mtag_b', 'mtag_d', 'mtag_h', 'mtag_s']
print("\t".join(output))
for i in range(0,22):
  mtag = mtag_list[i]
  pat = pat_list[i]
  migwas = migwas_list[i]
  hipo = hipo_list[i]
  
  mtag_sig = mtag[(mtag["mtag_pval1"] <= 5e-8) | (mtag["mtag_pval2"] <= 5e-8) | (mtag["mtag_pval3"] <= 5e-8) | (mtag["mtag_pval4"] <= 5e-8)].shape[0]
  mtag1 = mtag[mtag["mtag_pval1"] <= 5e-8].shape[0]
  mtag2 = mtag[mtag["mtag_pval2"] <= 5e-8].shape[0]
  mtag3 = mtag[mtag["mtag_pval3"] <= 5e-8].shape[0]
  mtag4 = mtag[mtag["mtag_pval4"] <= 5e-8].shape[0]

  pat_sig = pat[pat["PAT_Pvalue"] <= 5e-8].shape[0]
  pat1 = pat[(pat["PAT_Mvalue_1"] > 0.90) & (pat["PAT_Pvalue"] <= 5e-8)].shape[0]
  pat2 = pat[(pat["PAT_Mvalue_2"] > 0.90) & (pat["PAT_Pvalue"] <= 5e-8)].shape[0]
  pat3 = pat[(pat["PAT_Mvalue_3"] > 0.90) & (pat["PAT_Pvalue"] <= 5e-8)].shape[0]
  pat4 = pat[(pat["PAT_Mvalue_4"] > 0.90) & (pat["PAT_Pvalue"] <= 5e-8)].shape[0]

  pat12 = pat[(pat["PAT_Mvalue_1"] > 0.99) & (pat["PAT_Pvalue"] <= 5e-8)].shape[0]
  pat22 = pat[(pat["PAT_Mvalue_2"] > 0.99) & (pat["PAT_Pvalue"] <= 5e-8)].shape[0]
  pat32 = pat[(pat["PAT_Mvalue_3"] > 0.99) & (pat["PAT_Pvalue"] <= 5e-8)].shape[0]
  pat42 = pat[(pat["PAT_Mvalue_4"] > 0.99) & (pat["PAT_Pvalue"] <= 5e-8)].shape[0]

  migwas_sig = migwas[migwas["MIGWAS_Pvalue"] <= 5e-8].shape[0]
  migwas1 = migwas[(migwas["MIGWAS_Mvalue_1"] > 0.90) & (migwas["MIGWAS_Pvalue"] <= 5e-8)].shape[0]
  migwas2 = migwas[(migwas["MIGWAS_Mvalue_2"] > 0.90) & (migwas["MIGWAS_Pvalue"] <= 5e-8)].shape[0]
  migwas3 = migwas[(migwas["MIGWAS_Mvalue_3"] > 0.90) & (migwas["MIGWAS_Pvalue"] <= 5e-8)].shape[0]
  migwas4 = migwas[(migwas["MIGWAS_Mvalue_4"] > 0.90) & (migwas["MIGWAS_Pvalue"] <= 5e-8)].shape[0]

  migwas12 = migwas[(migwas["MIGWAS_Mvalue_1"] > 0.99) & (migwas["MIGWAS_Pvalue"] <= 5e-8)].shape[0]
  migwas22 = migwas[(migwas["MIGWAS_Mvalue_2"] > 0.99) & (migwas["MIGWAS_Pvalue"] <= 5e-8)].shape[0]
  migwas32 = migwas[(migwas["MIGWAS_Mvalue_3"] > 0.99) & (migwas["MIGWAS_Pvalue"] <= 5e-8)].shape[0]
  migwas42 = migwas[(migwas["MIGWAS_Mvalue_4"] > 0.99) & (migwas["MIGWAS_Pvalue"] <= 5e-8)].shape[0]

  hipo_sig = hipo[(hipo["pval.HIPOD1"] <= 5e-8) | (hipo["pval.HIPOD2"] <= 5e-8) | (hipo["pval.HIPOD3"] <= 5e-8) | (hipo["pval.HIPOD4"] <= 5e-8)].shape[0]
  hipo1 = hipo[(hipo["HIPO_Mvalue_0"] > 0.90) & ((hipo["pval.HIPOD1"] <= 5e-8) | (hipo["pval.HIPOD2"] <= 5e-8) | (hipo["pval.HIPOD3"] <= 5e-8) | (hipo["pval.HIPOD4"] <= 5e-8))].shape[0]
  hipo2 = hipo[(hipo["HIPO_Mvalue_1"] > 0.90) & ((hipo["pval.HIPOD1"] <= 5e-8) | (hipo["pval.HIPOD2"] <= 5e-8) | (hipo["pval.HIPOD3"] <= 5e-8) | (hipo["pval.HIPOD4"] <= 5e-8))].shape[0]
  hipo3 = hipo[(hipo["HIPO_Mvalue_2"] > 0.90) & ((hipo["pval.HIPOD1"] <= 5e-8) | (hipo["pval.HIPOD2"] <= 5e-8) | (hipo["pval.HIPOD3"] <= 5e-8) | (hipo["pval.HIPOD4"] <= 5e-8))].shape[0]
  hipo4 = hipo[(hipo["HIPO_Mvalue_3"] > 0.90) & ((hipo["pval.HIPOD1"] <= 5e-8) | (hipo["pval.HIPOD2"] <= 5e-8) | (hipo["pval.HIPOD3"] <= 5e-8) | (hipo["pval.HIPOD4"] <= 5e-8))].shape[0]

  hipo12 = hipo[(hipo["HIPO_Mvalue_0"] > 0.99) & ((hipo["pval.HIPOD1"] <= 5e-8) | (hipo["pval.HIPOD2"] <= 5e-8) | (hipo["pval.HIPOD3"] <= 5e-8) | (hipo["pval.HIPOD4"] <= 5e-8))].shape[0]
  hipo22 = hipo[(hipo["HIPO_Mvalue_1"] > 0.99) & ((hipo["pval.HIPOD1"] <= 5e-8) | (hipo["pval.HIPOD2"] <= 5e-8) | (hipo["pval.HIPOD3"] <= 5e-8) | (hipo["pval.HIPOD4"] <= 5e-8))].shape[0]
  hipo32 = hipo[(hipo["HIPO_Mvalue_2"] > 0.99) & ((hipo["pval.HIPOD1"] <= 5e-8) | (hipo["pval.HIPOD2"] <= 5e-8) | (hipo["pval.HIPOD3"] <= 5e-8) | (hipo["pval.HIPOD4"] <= 5e-8))].shape[0]
  hipo42 = hipo[(hipo["HIPO_Mvalue_3"] > 0.99) & ((hipo["pval.HIPOD1"] <= 5e-8) | (hipo["pval.HIPOD2"] <= 5e-8) | (hipo["pval.HIPOD3"] <= 5e-8) | (hipo["pval.HIPOD4"] <= 5e-8))].shape[0]

  output = [effect[i], causal_snps[i], pat_sig, pat1, pat2, pat3, pat4,
                                pat12, pat22, pat32, pat42,
                       migwas_sig, migwas1, migwas2, migwas3, migwas4,
                                migwas12, migwas22, migwas32, migwas42,
                       hipo_sig, hipo1, hipo2, hipo3, hipo4,
                                hipo12, hipo22, hipo32, hipo42,
                       mtag_sig, mtag1, mtag2, mtag3, mtag4]
  #print(output)
  output = [str(i) for i in output]
  print("\t".join(output))
