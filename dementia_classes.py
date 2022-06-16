from fn_cfg import *
import params as cfg

df = pd.read_excel (r'/Volumes/Backup Plus/EEG_Datasets/laurel_place/dataset/lp_clinical.xlsx')
data = (df[['ID','MMSE']]).to_numpy()

""""
def oneFourDigit(x):
    return '{:04d}'.format(x)

arr = []
for i in range(len(data)):
    arr.append(oneFourDigit(data[i,0]))
arr = np.array(arr)
data = np.vstack((arr, data[:,1])).T
"""

#   no dementia (ND)
idx_ND = np.where(data[:,1] >= 25)[0]
data_ND = data[:,0][idx_ND]

#   mild dementia (MILD)
idx_MILD = np.where(np.logical_and(data[:,1]>=20,data[:,1]<=24))
data_MILD = data[:,0][idx_MILD]

#   moderate dementia (MOD)
idx_MOD = np.where(np.logical_and(data[:,1]>=13,data[:,1]<=19))
data_MOD = data[:,0][idx_MOD]

#   severe dementia (SEVD)
idx_SEVD = np.where(data[:,1]<=12)
data_SEVD = data[:,0][idx_SEVD]

#   use the seperated IDs
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


