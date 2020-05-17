import pandas as pd
import sys
import numpy as np

direc=sys.argv[1]
pat_dir=sys.argv[2]
hipo_dir=sys.argv[3]
p = pd.read_csv(pat_dir, sep = ' ')
g = pd.read_csv(hipo_dir, sep = '\t')
shape = p.shape[0]

thresh = float(sys.argv[4])
p.loc[(p.PAT_Pvalue > thresh),"PAT_Mvalue_1"] = -1
p.loc[(p.PAT_Pvalue > thresh),"PAT_Mvalue_2"] = -1
p.loc[(p.PAT_Pvalue > thresh),"PAT_Mvalue_3"] = -1
p.loc[(p.PAT_Pvalue > thresh),"PAT_Mvalue_4"] = -1

g.loc[((g.pval_HIPOD1 > thresh) & (g.pval_HIPOD2 > thresh) & (g.pval_HIPOD3 > thresh) & (g.pval_HIPOD4 > thresh)),"HIPO_Mvalue_0"] = -1
g.loc[((g.pval_HIPOD1 > thresh) & (g.pval_HIPOD2 > thresh) & (g.pval_HIPOD3 > thresh) & (g.pval_HIPOD4 > thresh)),"HIPO_Mvalue_1"] = -1
g.loc[((g.pval_HIPOD1 > thresh) & (g.pval_HIPOD2 > thresh) & (g.pval_HIPOD3 > thresh) & (g.pval_HIPOD4 > thresh)),"HIPO_Mvalue_2"] = -1
g.loc[((g.pval_HIPOD1 > thresh) & (g.pval_HIPOD2 > thresh) & (g.pval_HIPOD3 > thresh) & (g.pval_HIPOD4 > thresh)),"HIPO_Mvalue_3"] = -1

m1 = pd.read_csv(direc+"/mtag_results_trait_1.txt", sep = '\t')
m2 = pd.read_csv(direc+"/mtag_results_trait_2.txt", sep = '\t')
m3 = pd.read_csv(direc+"/mtag_results_trait_3.txt", sep = '\t')
m4 = pd.read_csv(direc+"/mtag_results_trait_4.txt", sep = '\t')

mtag = m1["mtag_pval"].values
mtag = np.append(mtag,m2["mtag_pval"].values)
mtag = np.append(mtag,m3["mtag_pval"].values)
mtag = np.append(mtag,m4["mtag_pval"].values)

pat = p["PAT_Mvalue_1"].values
pat = np.append(pat,p["PAT_Mvalue_2"].values)
pat = np.append(pat,p["PAT_Mvalue_3"].values)
pat = np.append(pat,p["PAT_Mvalue_4"].values)

hipo = g["HIPO_Mvalue_0"].values
hipo = np.append(hipo,g["HIPO_Mvalue_1"].values)
hipo = np.append(hipo,g["HIPO_Mvalue_2"].values)
hipo = np.append(hipo,g["HIPO_Mvalue_3"].values)


#size = int((pat.shape[0] - 630000))
label = np.array([[1,1,1,1],[1,1,1,0],[1,0,1,0],[1,0,0,0],[0,0,1,0],[0,1,1,1],[0,1,0,1]])

total = np.tile([0,0,0,0],(int(630000),1))
shaper = 0.1*shape
shaper = shaper/7.0
for i in label:
  total = np.append(total,np.tile(i,(int(shaper),1)))
total = total.reshape((int(700000),4))
total = pd.DataFrame(total)

lab = total.iloc[:,0].values
lab = np.append(lab,total.iloc[:,1].values)
lab = np.append(lab,total.iloc[:,2].values)
lab = np.append(lab,total.iloc[:,3].values)

mtag_sort = pd.DataFrame({'truth':lab, 'MTAG':mtag})
pat_sort = pd.DataFrame({'truth':lab, 'PAT':pat})
hipo_sort = pd.DataFrame({'truth':lab, 'HIPO':hipo})

mtag_sort = mtag_sort.sort_values(by=['MTAG'])
mtag_sort = mtag_sort.reset_index(drop=True)
#print(mtag_sort.head())

pat_sort = pat_sort.sort_values(by=['PAT'], ascending = False)
pat_sort = pat_sort.reset_index(drop=True)
#print(pat_sort.head())

hipo_sort = hipo_sort.sort_values(by=['HIPO'], ascending = False)
hipo_sort = hipo_sort.reset_index(drop=True)
size=1000
pat_list = [pat_sort.loc[i:i+size-1,:] for i in range(0, pat_sort.shape[0],size)]
mtag_list = [mtag_sort.loc[i:i+size-1,:] for i in range(0, mtag_sort.shape[0],size)]
hipo_list = [hipo_sort.loc[i:i+size-1,:] for i in range(0,hipo_sort.shape[0],size)]

output = ['pat_sig', 'pat_mval1', 'pat_mval2', 'pat_mval3', 'pat_tp', 'pat_fp', 'mtag_sig', 'mtag_pval1', 'mtag_pval2', 'mtag_pval3', 'mtag_tp', 'mtag_fp', 'hipo_sig', 'hipo_mval1', 'hipo_mval2', 'hipo_mval3', 'hipo_tp', 'hipo_fp']
print("\t".join(output))
for i in range(0,len(pat_list)):
  mtag = mtag_list[i]
  pat = pat_list[i]
  hipo = hipo_list[i]

  pat_sig = pat[pat["PAT"] != -1].shape[0]
  pat = pat[pat["PAT"] != -1]
  if pat_sig > 0:
    pat_tp = pat[pat["truth"] == 1].shape[0]/pat_sig
    pat_fp = pat[pat["truth"] == 0].shape[0]/pat_sig
    pat_mval1 = pat.iloc[0]["PAT"]
    pat_mval2 = pat.iloc[-1]["PAT"]
    pat_mval3 = (pat_mval1 + pat_mval2)/2
  else:
    pat_tp = -1
    pat_fp = -1
    pat_mval1 = -1
    pat_mval2 = -1
    pat_mval3 = -1

  mtag_sig = mtag[mtag["MTAG"] <= thresh].shape[0]
  mtag = mtag[mtag["MTAG"] <= thresh]
  if mtag_sig > 0:
    mtag_tp = mtag[mtag["truth"] == 1].shape[0]/mtag_sig
    mtag_fp = mtag[mtag["truth"] == 0].shape[0]/mtag_sig
    mtag_pval1 = mtag.iloc[0]["MTAG"]
    mtag_pval2 = mtag.iloc[-1]["MTAG"]
    mtag_pval3 = (mtag_pval1 + mtag_pval2)/2
  else:
    mtag_tp = -1
    mtag_fp = -1
    mtag_pval1 = -1
    mtag_pval2 = -1
    mtag_pval3 = -1

  hipo_sig = hipo[hipo["HIPO"] != -1].shape[0]
  hipo = hipo[hipo["HIPO"] != -1]
  if hipo_sig > 0: 
    hipo_tp = hipo[hipo["truth"] == 1].shape[0]/hipo_sig
    hipo_fp = hipo[hipo["truth"] == 0].shape[0]/hipo_sig
    hipo_mval1 = hipo.iloc[0]["HIPO"]
    hipo_mval2 = hipo.iloc[-1]["HIPO"]
    hipo_mval3 = (hipo_mval1 + hipo_mval2)/2
  else:
    hipo_tp = -1
    hipo_fp = -1
    hipo_mval1 = -1
    hipo_mval2 = -1
    hipo_mval3 = -1

  output = [pat_sig, pat_mval1, pat_mval2, pat_mval3, pat_tp, pat_fp, mtag_sig, mtag_pval1, mtag_pval2, mtag_pval3, mtag_tp, mtag_fp, hipo_sig, hipo_mval1, hipo_mval2, hipo_mval3, hipo_tp, hipo_fp]
  output = [str(i) for i in output]
  print("\t".join(output))
