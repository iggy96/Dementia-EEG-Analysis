"""
Joshua Ighalo 12/06/2022

dementia study scans sorter

this script sorts the scans into the following categories:
1. severe dementia
2. moderate dementia
3. mild dementia
4. no dementia

Input: .xlsx file with the following columns: scan ID and MMSE
       



"""
from fn_cfg import *
import params as cfg

df_1 = pd.read_excel (r'/Users/joshuaighalo/Downloads/EEG_Datasets/laurel_place/dataset/lp_clinical.xlsx')
df_2 = pd.read_excel(r'/Users/joshuaighalo/Downloads/EEG_Datasets/laurel_place/dataset/Scan_details.xlsx')
data_1 = (df_1[['ID','MMSE']]).to_numpy()
data_2 = (df_2[['Baseline','4 months','8 months']]).to_numpy()

#   no dementia (ND)
idx_ND = np.where(data_1[:,1] >= 25)[0]
data_ND = data_1[:,0][idx_ND]

#   mild dementia (MILD)
idx_MILD = np.where(np.logical_and(data_1[:,1]>=20,data_1[:,1]<=24))
data_MILD = data_1[:,0][idx_MILD]

#   moderate dementia (MOD)
idx_MOD = np.where(np.logical_and(data_1[:,1]>=13,data_1[:,1]<=19))
data_MOD = data_1[:,0][idx_MOD]

#   severe dementia (SEVD)
idx_SEVD = np.where(data_1[:,1]<=12)
data_SEVD = data_1[:,0][idx_SEVD]

#   seperate scan IDs for each dementia type
def oneFourDigit(x):
    return '{:04d}'.format(x)

arr = []
for i in range(len(data_ND)):
    arr.append(oneFourDigit(data_ND[i]))
data_ND = np.array(arr)
arr = []
for i in range(len(data_MILD)):
    arr.append(oneFourDigit(data_MILD[i]))
data_MILD = np.array(arr)
arr = []
for i in range(len(data_MOD)):
    arr.append(oneFourDigit(data_MOD[i]))
data_MOD = np.array(arr)
arr = []
for i in range(len(data_SEVD)):
    arr.append(oneFourDigit(data_SEVD[i]))
data_SEVD = np.array(arr)

#   import the names of all the scans folders
import os
search_path = r'/Users/joshuaighalo/Downloads/EEG_Datasets/laurel_place/dataset'
root, dirs, files = next(os.walk(search_path), ([],[],[]))
scanFolders = np.array(dirs)

#   match scan IDs from csv with scan folders
scanFolders_Init = np.array([w[0:4] for w in scanFolders])
idx_ND = np.where(np.isin(scanFolders_Init, data_ND))[0]
scans_no_dementia = scanFolders[idx_ND]
idx_MILD = np.where(np.isin(scanFolders_Init, data_MILD))[0]
scans_mild_dementia = scanFolders[idx_MILD]
idx_MOD = np.where(np.isin(scanFolders_Init, data_MOD))[0]
scans_mod_dementia = scanFolders[idx_MOD]
idx_SEVD = np.where(np.isin(scanFolders_Init, data_SEVD))[0]
scans_sev_dementia = scanFolders[idx_SEVD]

#   use dates from scan details to sort scans into timepoints 
#   transform date from yy-mm-dd to ddmmyyyy to match filing system of scans
data_2[np.isnan(data_2)] = 0    # replace nan with unused date

def ext(data):
    char = datetime.strptime(str(data)[0:10], '%Y-%m-%d').strftime('%d%m%Y')
    return char

tp_baseline = []
tp_4months = []
tp_8months = []
for i in range(len(data_2)):
    tp_baseline.append(ext(data_2[i,0]))
    tp_4months.append(ext(data_2[i,1]))
    tp_8months.append(ext(data_2[i,2]))

#   use the timepoints to sort scans per timepoints
def remove_cruft(s):
    return s[7:-5]

tp_ND = np.array([remove_cruft(s) for s in scans_no_dementia])
tp_MOD = np.array([remove_cruft(s) for s in scans_mod_dementia])
tp_MILD = np.array([remove_cruft(s) for s in scans_mild_dementia])
tp_SEV = np.array([remove_cruft(s) for s in scans_sev_dementia])


#   extract baseline, 4-months and 8-months for the no dementia scans
idx_B_ND = 
