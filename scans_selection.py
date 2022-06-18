from fn_cfg import *
import params as cfg
import dementia_classes


#   baseline scans for different dementia classes
base_ND = dementia_classes.scansND_B
base_MILD = dementia_classes.scansMILD_B
base_MOD = dementia_classes.scansMOD_B
base_SD = dementia_classes.scansSEVD_B


#   4-months scans for different dementia classes
four_ND = dementia_classes.scansND_4
four_MILD = dementia_classes.scansMILD_4
four_MOD = dementia_classes.scansMOD_4
four_SD = dementia_classes.scansSEVD_4
four_run1_SD = dementia_classes.run1_scansSEVD_4
four_run2_SD = dementia_classes.run2_scansSEVD_4

#   8-months scans for different dementia classes
eight_ND = dementia_classes.scansND_8
eight_MILD = dementia_classes.scansMILD_8
eight_MOD = dementia_classes.scansMOD_8
eight_SD = dementia_classes.scansSEVD_8


def eegSignalQuality(device_version,scan_ID,threshold,channel_name,local_path,dispIMG):
    """
    """
    chans = dict(Fz=0,Cz=1,Pz=2)
    version,filename,localPath = device_version,scan_ID,local_path 
    device = importFile.neurocatch()
    fileObjects = device.init(version,filename,localPath,dispIMG=False)
    rawEEG = fileObjects[0]
    rawEOG = fileObjects[1]
    filtering = filters()
    adaptiveFilterOutput = filtering.adaptive(rawEEG,rawEOG)
    notchFilterOutput = filtering.notch(adaptiveFilterOutput,line,fs)
    bandPassFilterOutput = filtering.butterBandPass(notchFilterOutput,lowcut=0.1,highcut=5,fs=cfg.fs)
    chan_idx = chans[channel_name]
    pkScore = quality_p2p(bandPassFilterOutput)
    pkScore = pkScore[chan_idx]
    if pkScore <= threshold:
        if dispIMG == True:
            print(scan_ID,"quality: good")
        else:
            pass
        score = pkScore
        scan = scan_ID
    elif pkScore > threshold:
        if dispIMG == True:
            print(scan_ID,"quality: bad")
        else:
            pass
        score = pkScore
        scan = 'rejected'
    if math.isnan(pkScore)==True:
        if dispIMG == True:
            print(scan_ID,"quality: bad")
        else:
            pass
        score = pkScore
        scan = 'rejected'
    return scan,score



threshold = 1000
localPath = '/Users/joshuaighalo/Downloads/EEG_Datasets/laurel_place/dataset'
version = 1.0
channel = 'Cz'
dispIMG = True

run_IDs_1 = four_run1_SD
run_IDs_2 = four_run2_SD

status_runs_1 = []
status_runs_2 = []
scores_runs_1 = []
scores_runs_2 = []
print("Number of scans in run 1: ",len(run_IDs_1))
print("Number of scans in run 2: ",len(run_IDs_2))
for i in range(len(run_IDs_1)):
    scan_ID = run_IDs_1[i]
    scan,score = eegSignalQuality(version,scan_ID,threshold,channel,localPath,dispIMG)
    status_runs_1.append(scan)
    scores_runs_1.append(score)

for i in range(len(run_IDs_2)):
    scan_ID = run_IDs_2[i]
    scan,score = eegSignalQuality(version,scan_ID,threshold,channel,localPath,dispIMG)
    status_runs_2.append(scan)
    scores_runs_2.append(score)

idx_rejected_1 = [i for i,val in enumerate(status_runs_1) if val=='rejected']
idx_accepted_1 = [i for i,val in enumerate(status_runs_1) if val!='rejected']
idx_rejected_2 = [i for i,val in enumerate(status_runs_2) if val=='rejected']
idx_accepted_2 = [i for i,val in enumerate(status_runs_2) if val!='rejected']
runIDs_1_Rej = [run_IDs_1[i] for i in idx_rejected_1]
runIDs_1_Acc = [run_IDs_1[i] for i in idx_accepted_1]
runIDs_2_Rej = [run_IDs_2[i] for i in idx_rejected_2]
runIDs_2_Acc = [run_IDs_2[i] for i in idx_accepted_2]
scores_1_Acc = [scores_runs_1[i] for i in idx_accepted_1]
scores_1_Rej = [scores_runs_1[i] for i in idx_rejected_1]
scores_2_Acc = [scores_runs_2[i] for i in idx_accepted_2]
scores_2_Rej = [scores_runs_2[i] for i in idx_rejected_2]
print('\n')
print("Number of accepted scans in run 1: ",len(runIDs_1_Acc))
print("Number of rejected scans in run 1: ",len(runIDs_1_Rej))
print('\n')
print("Number of accepted scans in run 2: ",len(runIDs_2_Acc))
print("Number of rejected scans in run 2: ",len(runIDs_2_Rej))

#char_runIDs_Acc_1 = [x[:-16] for x in runIDs_1_Acc]
#char_runIDs_Acc_2 = [x[:-16] for x in runIDs_2_Acc]

merge_acc_runIDs = runIDs_1_Acc + runIDs_2_Acc
merge_acc_scores = scores_1_Acc + scores_2_Acc

char_merge_acc_runIDs = [x[:-16] for x in merge_acc_runIDs]
idx_non_dups = [idx_ for idx_, item in enumerate(char_merge_acc_runIDs) if item not in char_merge_acc_runIDs[:idx_]]


