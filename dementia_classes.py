from fn_cfg import *
import params as cfg

df = pd.read_excel (r'/Volumes/Backup Plus/EEG_Datasets/laurel_place/dataset/lp_clinical.xlsx')
df_new = (df[['ID','MMSE']]).to_numpy()

#   no dementia
df_no_dementia = df_new[df_new[:,1] >= 25]
#   mild dementia
df_mild_dementia = df_new[df_new[:,1] 
