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
"""
df_1 = pd.read_excel (r'/Users/joshuaighalo/Downloads/EEG_Datasets/laurel_place/dataset/lp_clinical.xlsx')
df_2 = pd.read_excel(r'/Users/joshuaighalo/Downloads/EEG_Datasets/laurel_place/dataset/Scan_details.xlsx')
data_1 = (df_1[['ID','MMSE']]).to_numpy()
data_2 = (df_2[['Baseline','4 months','8 months']]).to_numpy()

#   remove scans of participants who died
died = [152,280,409]
idx_died = np.where(np.isin(data_1[:,0], died))[0]
data_1 = np.delete(data_1,idx_died,axis=0)
data_2 = np.delete(data_2,idx_died,axis=0)

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
scans_ND = scanFolders[idx_ND]
idx_MILD = np.where(np.isin(scanFolders_Init, data_MILD))[0]
scans_MILD = scanFolders[idx_MILD]
idx_MOD = np.where(np.isin(scanFolders_Init, data_MOD))[0]
scans_MOD = scanFolders[idx_MOD]
idx_SEVD = np.where(np.isin(scanFolders_Init, data_SEVD))[0]
scans_SEVD = scanFolders[idx_SEVD]

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

#   remove nan cells (in the form of 01011970) from timepoint arrays
idx_nan_B = np.where(np.isin(tp_baseline, "01011970"))[0]
tp_baseline = list(np.delete(tp_baseline,idx_nan_B,axis=0))
idx_nan_4 = np.where(np.isin(tp_4months, "01011970"))[0]
tp_4months = list(np.delete(tp_4months,idx_nan_4,axis=0))
idx_nan_8 = np.where(np.isin(tp_8months, "01011970"))[0]
tp_8months = list(np.delete(tp_8months,idx_nan_8,axis=0))

#   use the timepoints to sort scans per timepoints
def remove_cruft(s):
    return s[7:-5]

tp_ND = np.array([remove_cruft(s) for s in scans_ND])
tp_MOD = np.array([remove_cruft(s) for s in scans_MOD])
tp_MILD = np.array([remove_cruft(s) for s in scans_MILD])
tp_SEV = np.array([remove_cruft(s) for s in scans_SEVD])


#   extract baseline, 4-months and 8-months for the no dementia scans
idx_ND_B = np.where(np.isin(tp_ND, np.array(tp_baseline)))[0]
idx_ND_4 = np.where(np.isin(tp_ND, tp_4months))[0]
idx_ND_8 = np.where(np.isin(tp_ND, tp_8months))[0]
scans_ND_B = scans_ND[idx_ND_B]
scans_ND_4 = scans_ND[idx_ND_4]
scans_ND_8 = scans_ND[idx_ND_8]

#   extract baseline, 4-months and 8-months for the mild dementia scans
idx_MILD_B = np.where(np.isin(tp_MILD, tp_baseline))[0]
idx_MILD_4 = np.where(np.isin(tp_MILD, tp_4months))[0]
idx_MILD_8 = np.where(np.isin(tp_MILD, tp_8months))[0]
scans_MILD_B = scans_MILD[idx_MILD_B]
scans_MILD_4 = scans_MILD[idx_MILD_4]
scans_MILD_8 = scans_MILD[idx_MILD_8]

#   extract baseline, 4-months and 8-months for the moderate dementia scans
idx_MOD_B = np.where(np.isin(tp_MOD, tp_baseline))[0]
idx_MOD_4 = np.where(np.isin(tp_MOD, tp_4months))[0]
idx_MOD_8 = np.where(np.isin(tp_MOD, tp_8months))[0]
scans_MOD_B = scans_MOD[idx_MOD_B]
scans_MOD_4 = scans_MOD[idx_MOD_4]
scans_MOD_8 = scans_MOD[idx_MOD_8]

#   extract baseline, 4-months and 8-months for the severe dementia scans
idx_SEVD_B = np.where(np.isin(tp_SEV, tp_baseline))[0]
idx_SEVD_4 = np.where(np.isin(tp_SEV, tp_4months))[0]
idx_SEVD_8 = np.where(np.isin(tp_SEV, tp_8months))[0]
scans_SEVD_B = scans_SEVD[idx_SEVD_B]
scans_SEVD_4 = scans_SEVD[idx_SEVD_4]
scans_SEVD_8 = scans_SEVD[idx_SEVD_8]
"""






df_1 = pd.read_excel (r'/Users/joshuaighalo/Downloads/EEG_Datasets/laurel_place/dataset/lp_clinical.xlsx')
df_2 = pd.read_excel(r'/Users/joshuaighalo/Downloads/EEG_Datasets/laurel_place/dataset/Scan_details.xlsx')
data_1 = (df_1[['ID','MMSE']]).to_numpy()
data_2 = (df_2[['Baseline','4 months','8 months']]).to_numpy()

#   remove scans of participants who died
died = [152,280,409]
idx_died = np.where(np.isin(data_1[:,0], died))[0]
data_1a = np.delete(data_1,idx_died,axis=0)
data_1 = np.delete(data_1,idx_died,axis=0)
data_2 = np.delete(data_2,idx_died,axis=0)


#   seperate scan IDs for each dementia type
def oneFourDigit(x):
    return '{:04d}'.format(x)

arr = []
for i in range(len(data_1)):
    arr.append(oneFourDigit(data_1[:,0][i]))
data_10 = np.array(arr)
data_1 = np.vstack((data_10,data_1[:,1])).T

#   create an array form of the timepoints csv
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
tp_baseline = np.array(tp_baseline)
tp_4months = np.array(tp_4months)
tp_8months = np.array(tp_8months)

#   merge data_1 and timepoints
def stringMerge(x,y):
    return x + y
dataBase = []
data4months = []
data8months = []
for i in range(len(data_1)):
    dataBase.append(stringMerge(data_1[:,0][i], tp_baseline[i]))
    data4months.append(stringMerge(data_1[:,0][i], tp_4months[i]))
    data8months.append(stringMerge(data_1[:,0][i], tp_8months[i]))

#   import dataset scan from directory
search_path = r'/Users/joshuaighalo/Downloads/EEG_Datasets/laurel_place/dataset'
root, dirs, files = next(os.walk(search_path), ([],[],[]))
scans = np.array(dirs)

#   remove "_runNumber_" from scan names and remove last five 
#   characters from scan names to get ID
def kill_char(string, indices): 
    #   e.g., data, indexes = "Welcome", {1, 3, 5}
    data, indexes = string,indices
    new_str = "".join([char for idx, char in enumerate(data) if idx not in indexes])
    new_str = new_str[:-5]
    return new_str

dataScans = []
for i in range(len(scans)):
    dataScans.append(kill_char(scans[i],{4,5,6}))
dataScans = np.array(dataScans)

#   baseline scans in the dataset and their positions in the dataset
#   Note: index positions in dataScans = index positions in scans
idx_B = np.where(np.isin(dataScans, dataBase))[0]
dataScans_B = dataScans[idx_B]
idx_4 = np.where(np.isin(dataScans, data4months))[0]
dataScans_4 = dataScans[idx_4]
idx_8 = np.where(np.isin(dataScans, data8months))[0]
dataScans_8 = dataScans[idx_8]

"""
-   Apply mmse threshold as defined by alz.org to classify participants into dementia classes
-   len and indices of data_1a == len and indices of data_1
-   no dementia (ND)
"""
#   no dementia (ND)
idx_ND = np.where(data_1a[:,1] >= 25)[0]
data_ND = data_1[:,0][idx_ND]

#   mild dementia (MILD)
idx_MILD = np.where(np.logical_and(data_1a[:,1]>=20,data_1a[:,1]<=24))
data_MILD = data_1[:,0][idx_MILD]

#   moderate dementia (MOD)
idx_MOD = np.where(np.logical_and(data_1a[:,1]>=13,data_1a[:,1]<=19))
data_MOD = data_1[:,0][idx_MOD]

#   severe dementia (SEVD)
idx_SEVD = np.where(data_1a[:,1]<=12)
data_SEVD = data_1[:,0][idx_SEVD]

"""
find the scans of the classes of dementia within the timepoints scans
"""

#   scans of no dementia across timepoints
def test(data_class, dataScans_TP,dataScans,scans):
    idx_class_TP = np.where(np.isin(([x[:-8] for x in dataScans_TP]), data_class))[0]
    init_sname = dataScans_TP[idx_class_TP]
    idx_scansClass_TP = np.where(np.isin(dataScans, init_sname))[0]
    scansClass_TP = scans[idx_scansClass_TP]  
    return scansClass_TP 

# no dementia (ND) scans at baseline
scansND_B = []
for i in range(len(data_ND)):
    scansND_B.append(test(data_ND[i], dataScans_B,dataScans,scans))
#print("scansND_B: ", scansND_B)
#idx_ND_B = np.where(np.isin(([x[:-8] for x in dataScans_B]), data_ND))[0]
#init_sname = dataScans_B[idx_ND_B]
#idx_scansND_B = np.where(np.isin(dataScans, init_sname))[0]
#scansND_B = scans[idx_scansND_B]


# no dementia (ND) scans at 4 months
scansND_4 = []
for i in range(len(data_ND)):
    scansND_4.append(test(data_ND[i], dataScans_4,dataScans,scans))
#print("scansND_4: ", scansND_4)
#idx_ND_4 = np.where(np.isin(([x[:-8] for x in dataScans_4]), data_ND))[0]
#init_sname = dataScans_4[idx_ND_4]
#idx_scansND_4 = np.where(np.isin(dataScans, init_sname))[0]
#scansND_4 = scans[idx_scansND_4]


# no dementia (ND) scans at 8 months
scansND_8 = []
for i in range(len(data_ND)):
    scansND_8.append(test(data_ND[i], dataScans_8,dataScans,scans))
#print("scansND_8: ", scansND_8)
#idx_ND_8 = np.where(np.isin(([x[:-8] for x in dataScans_8]), data_ND))[0]
#init_sname = dataScans_8[idx_ND_8]
#idx_scansND_8 = np.where(np.isin(dataScans, init_sname))[0]
#scansND_8 = scans[idx_scansND_8]


#  scans of mild dementia across timepoints

# mild dementia (MILD) scans at baseline
scansMILD_B = []
for i in range(len(data_MILD)):
    scansMILD_B.append(test(data_MILD[i], dataScans_B,dataScans,scans))

#idx_MILD_B = np.where(np.isin(([x[:-8] for x in dataScans_B]), data_MILD))[0]
#init_sname = dataScans_B[idx_MILD_B]
#idx_scansMILD_B = np.where(np.isin(dataScans, init_sname))[0]
#scansMILD_B = scans[idx_scansMILD_B]

# mild dementia (MILD) scans at 4 months
scansMILD_4 = []
for i in range(len(data_MILD)):
    scansMILD_4.append(test(data_MILD[i], dataScans_4,dataScans,scans))

#idx_MILD_4 = np.where(np.isin(([x[:-8] for x in dataScans_4]), data_MILD))[0]
#init_sname = dataScans_4[idx_MILD_4]
#idx_scansMILD_4 = np.where(np.isin(dataScans, init_sname))[0]
#scansMILD_4 = scans[idx_scansMILD_4]

# mild dementia (MILD) scans at 8 months
scansMILD_8 = []
for i in range(len(data_MILD)):
    scansMILD_8.append(test(data_MILD[i], dataScans_8,dataScans,scans))

#idx_MILD_8 = np.where(np.isin(([x[:-8] for x in dataScans_8]), data_MILD))[0]
#init_sname = dataScans_8[idx_MILD_8]
#idx_scansMILD_8 = np.where(np.isin(dataScans, init_sname))[0]
#scansMILD_8 = scans[idx_scansMILD_8]

#  scans of moderate dementia across timepoints

# moderate dementia (MOD) scans at baseline
scansMOD_B = []
for i in range(len(data_MOD)):
    scansMOD_B.append(test(data_MOD[i], dataScans_B,dataScans,scans))

#idx_MOD_B = np.where(np.isin(([x[:-8] for x in dataScans_B]), data_MOD))[0]
#init_sname = dataScans_B[idx_MOD_B]
#idx_scansMOD_B = np.where(np.isin(dataScans, init_sname))[0]
#scansMOD_B = scans[idx_scansMOD_B]

# moderate dementia (MOD) scans at 4 months
scansMOD_4 = []
for i in range(len(data_MOD)):
    scansMOD_4.append(test(data_MOD[i], dataScans_4,dataScans,scans))

#idx_MOD_4 = np.where(np.isin(([x[:-8] for x in dataScans_4]), data_MOD))[0]
#init_sname = dataScans_4[idx_MOD_4]
#idx_scansMOD_4 = np.where(np.isin(dataScans, init_sname))[0]
#scansMOD_4 = scans[idx_scansMOD_4]

# moderate dementia (MOD) scans at 8 months
scansMOD_8 = []
for i in range(len(data_MOD)):
    scansMOD_8.append(test(data_MOD[i], dataScans_8,dataScans,scans))

#idx_MOD_8 = np.where(np.isin(([x[:-8] for x in dataScans_8]), data_MOD))[0]
#init_sname = dataScans_8[idx_MOD_8]
#idx_scansMOD_8 = np.where(np.isin(dataScans, init_sname))[0]
#scansMOD_8 = scans[idx_scansMOD_8]

#  scans of severe dementia across timepoints

# severe dementia (SEVD) scans at baseline
scansSEVD_B = []
for i in range(len(data_SEVD)):
    scansSEVD_B.append(test(data_SEVD[i], dataScans_B,dataScans,scans))

#idx_SEVD_B = np.where(np.isin(([x[:-8] for x in dataScans_B]), data_SEVD))[0]
#init_sname = dataScans_B[idx_SEVD_B]
#idx_scansSEVD_B = np.where(np.isin(dataScans, init_sname))[0]
#scansSEVD_B = scans[idx_scansSEVD_B]

# severe dementia (SEVD) scans at 4 months
scansSEVD_4 = []
for i in range(len(data_SEVD)):
    scansSEVD_4.append(test(data_SEVD[i], dataScans_4,dataScans,scans))

#idx_SEVD_4 = np.where(np.isin(([x[:-8] for x in dataScans_4]), data_SEVD))[0]
#init_sname = dataScans_4[idx_SEVD_4]
#idx_scansSEVD_4 = np.where(np.isin(dataScans, init_sname))[0]
#scansSEVD_4 = scans[idx_scansSEVD_4]

# severe dementia (SEVD) scans at 8 months
scansSEVD_8 = []
for i in range(len(data_SEVD)):
    scansSEVD_8.append(test(data_SEVD[i], dataScans_8,dataScans,scans))
    
#idx_SEVD_8 = np.where(np.isin(([x[:-8] for x in dataScans_8]), data_SEVD))[0]
#init_sname = dataScans_8[idx_SEVD_8]
#idx_scansSEVD_8 = np.where(np.isin(dataScans, init_sname))[0]
#scansSEVD_8 = scans[idx_scansSEVD_8]









