"""
Input:  scans (variables) classfified by dementia class and timepoint from dementia_classes.py
Output: scans with the best eeg quality to represent each participants within a class

Note: -     each participant within the dementia classes have mutiple runs taken per timepoint 
            i.e., some participants have eeg taken from them twice per timepoint while it varies for others
      -     code is designed based on two runs and not three timepoints
"""
from fn_cfg import *
import params as cfg
import dementia_classes




print("NOTES:")
print("Pre-Quality Processing: Participants of each group have multiple scans per timepoints")
print("Post-Quality Processing: Participants of each group are allocated just one scan")


def scanSelection(device_version,run_IDs_1,run_IDs_2,threshold,channel_name,local_path,dem_class,timepoint,dispSubStats,dispMainStats):

    device_version,run_IDs_1,run_IDs_2,threshold,channel_name,local_path,dem_class,timepoint,dispSubStats,dispMainStats = device_version,run_IDs_1,run_IDs_2,threshold,channel_name,local_path,dem_class,timepoint,dispSubStats,dispMainStats
    
    
    def eegSignalQuality(device_version,scan_ID,threshold,channel_name,local_path,dispSubStats):
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
            if dispSubStats == True:
                print(scan_ID,"quality: good")
            else:
                pass
            score = pkScore
            scan = scan_ID
        elif pkScore > threshold:
            if dispSubStats == True:
                print(scan_ID,"quality: bad")
            else:
                pass
            score = pkScore
            scan = 'rejected'
        if math.isnan(pkScore)==True:
            if dispSubStats == True:
                print(scan_ID,"quality: bad")
            else:
                pass
            score = pkScore
            scan = 'rejected'
        return scan,score


    print("Processing for",dem_class,'at',timepoint,'.............................................................')
    char_run_IDs_1 = [x[:-16] for x in run_IDs_1]
    char_run_IDs_2 = [x[:-16] for x in run_IDs_2]
    char_runs = char_run_IDs_1 + char_run_IDs_2
    no_dup_chars = list(dict.fromkeys(char_runs))

    if dispMainStats == True:
        print("Number of",dem_class,"participants pre-quality assessment:",len(no_dup_chars))
    else:
        pass

    if dispSubStats==True:
        print("Number of scans in run 1: ",len(run_IDs_1))
        print("Number of scans in run 2: ",len(run_IDs_2))
    else:
        pass

    status_runs_1 = []
    status_runs_2 = []
    scores_runs_1 = []
    scores_runs_2 = []
    for i in range(len(run_IDs_1)):
        scan_ID = run_IDs_1[i]
        scan,score = eegSignalQuality(version,scan_ID,threshold,channel,localPath,dispSubStats)
        status_runs_1.append(scan)
        scores_runs_1.append(score)

    for i in range(len(run_IDs_2)):
        scan_ID = run_IDs_2[i]
        scan,score = eegSignalQuality(version,scan_ID,threshold,channel,localPath,dispSubStats)
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

    if dispSubStats==True:
        print('\n')
        print("Number of accepted scans in run 1: ",len(runIDs_1_Acc))
        print("Number of rejected scans in run 1: ",len(runIDs_1_Rej))
        print('\n')
        print("Number of accepted scans in run 2: ",len(runIDs_2_Acc))
        print("Number of rejected scans in run 2: ",len(runIDs_2_Rej))
    else:
        pass

    merge_acc_runIDs = runIDs_1_Acc + runIDs_2_Acc
    merge_acc_scores = scores_1_Acc + scores_2_Acc

    char_merge_acc_runIDs = [x[:-16] for x in merge_acc_runIDs]

    none_duplicate_charIDs = [unique for unique in (list(set(char_merge_acc_runIDs))) if char_merge_acc_runIDs.count(unique) == 1]
    d = {item: idx for idx, item in enumerate(char_merge_acc_runIDs)}
    items_to_find = none_duplicate_charIDs
    idx_none_duplicate_charIDs =[d.get(item) for item in items_to_find]
    none_duplicate_runIDs = [merge_acc_runIDs[i] for i in idx_none_duplicate_charIDs]

    new_acc_runIDs = [v for i, v in enumerate(merge_acc_runIDs) if i not in idx_none_duplicate_charIDs]
    new_acc_scores = [v for i, v in enumerate(merge_acc_scores) if i not in idx_none_duplicate_charIDs]

    rolling_scores = slidingWindow(np.array(new_acc_scores),int(len(new_acc_runIDs)/2),int(len(new_acc_runIDs)/2)).T
    rolling_runIDs = slidingWindow(np.array(new_acc_runIDs),int(len(new_acc_runIDs)/2),int(len(new_acc_runIDs)/2)).T
    idx_bestPeaks = np.argmin(rolling_scores,axis=1)
    bestPeaks_runIDs = [rolling_runIDs[i][j] for i,j in zip(range(len(rolling_runIDs)),idx_bestPeaks)]
    accepted_runIDs = bestPeaks_runIDs+none_duplicate_runIDs

    if dispMainStats == True:
        print("Number of",dem_class,"participants post-quality assessment:",len(accepted_runIDs))
    else:
        pass
    print('\n')
    return accepted_runIDs




#   baseline scans for different dementia classes
base_run1_ND = dementia_classes.run1_scansND_B
base_run2_ND = dementia_classes.run2_scansND_B
base_run1_MILD = dementia_classes.run1_scansMILD_B
base_run2_MILD = dementia_classes.run2_scansMILD_B
base_run1_MOD = dementia_classes.run1_scansMOD_B
base_run2_MOD = dementia_classes.run2_scansMOD_B
base_run1_SEVD = dementia_classes.run1_scansSEVD_B
base_run2_SEVD = dementia_classes.run2_scansSEVD_B

#   4-months scans for different dementia classes
four_run1_ND = dementia_classes.run1_scansND_4
four_run2_ND = dementia_classes.run2_scansND_4
four_run1_MILD = dementia_classes.run1_scansMILD_4
four_run2_MILD = dementia_classes.run2_scansMILD_4
four_run1_MOD = dementia_classes.run1_scansMOD_4
four_run2_MOD = dementia_classes.run2_scansMOD_4
four_run1_SEVD = dementia_classes.run1_scansSEVD_4
four_run2_SEVD = dementia_classes.run2_scansSEVD_4

#   8-months scans for different dementia classes
eight_run1_ND = dementia_classes.run1_scansND_8
eight_run2_ND = dementia_classes.run2_scansND_8
eight_run1_MILD = dementia_classes.run1_scansMILD_8
eight_run2_MILD = dementia_classes.run2_scansMILD_8
eight_run1_MOD = dementia_classes.run1_scansMOD_8
eight_run2_MOD = dementia_classes.run2_scansMOD_8
eight_run1_SEVD = dementia_classes.run1_scansSEVD_8
eight_run2_SEVD = dementia_classes.run2_scansSEVD_8


threshold = 300
localPath = '/Users/joshuaighalo/Downloads/EEG_Datasets/laurel_place/dataset'
version = 1.0
channel = 'Cz'
dispSubStats = False
dispMainStats = True

baseND = scanSelection(version,base_run1_ND,base_run2_ND,threshold,channel,localPath,
                        'no dementia','baseline',dispSubStats,dispMainStats)
baseMILD = scanSelection(version,base_run1_MILD,base_run2_MILD,threshold,channel,localPath,
                        'mild dementia','baseline',dispSubStats,dispMainStats)
baseMOD = scanSelection(version,base_run1_MOD,base_run2_MOD,threshold,channel,localPath,
                        'moderate dementia','baseline',dispSubStats,dispMainStats)
baseSEVD = scanSelection(version,base_run1_SEVD,base_run2_SEVD,threshold,channel,localPath,
                        'severe dementia','baseline',dispSubStats,dispMainStats)
fourND = scanSelection(version,four_run1_ND,four_run2_ND,threshold,channel,localPath,
                        'no dementia','4-months',dispSubStats,dispMainStats)
fourMILD = scanSelection(version,four_run1_MILD,four_run2_MILD,threshold,channel,localPath,
                        'mild dementia','4-months',dispSubStats,dispMainStats)
fourMOD = scanSelection(version,four_run1_MOD,four_run2_MOD,threshold,channel,localPath,
                        'moderate dementia','4-months',dispSubStats,dispMainStats)
fourSEVD = scanSelection(version,four_run1_SEVD,four_run2_SEVD,threshold,channel,localPath,
                        'severe dementia','4-months',dispSubStats,dispMainStats)
eightND = scanSelection(version,eight_run1_ND,eight_run2_ND,threshold,channel,localPath,
                        'no dementia','8-months',dispSubStats,dispMainStats)
eightMILD = scanSelection(version,eight_run1_MILD,eight_run2_MILD,threshold,channel,localPath,
                        'mild dementia','8-months',dispSubStats,dispMainStats)
eightMOD = scanSelection(version,eight_run1_MOD,eight_run2_MOD,threshold,channel,localPath,
                        'moderate dementia','8-months',dispSubStats,dispMainStats)
eightSEVD = scanSelection(version,eight_run1_SEVD,eight_run2_SEVD,threshold,channel,localPath,
                        'severe dementia','8-months',dispSubStats,dispMainStats)
