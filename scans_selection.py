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

#   8-months scans for different dementia classes
eight_ND = dementia_classes.scansND_8
eight_MILD = dementia_classes.scansMILD_8
eight_MOD = dementia_classes.scansMOD_8
eight_SD = dementia_classes.scansSEVD_8


def eegSignalQuality(device_version,scan_ID,threshold,channel_name,local_path):
    """
    """
    print(scan_ID)
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
        scan = scan_ID
    elif pkScore > threshold:
        scan = 'rejected'
    if math.isnan(pkScore)==True:
        scan = 'rejected'
    return scan



threshold = 150
localPath = '/Users/joshuaighalo/Downloads/EEG_Datasets/laurel_place/dataset'
version = 1.0
channel = 'Fz'

run_IDs_1 = four_run1_SD

score_runs = []
for i in range(len(run_IDs_1)):
    score_runs.append(eegSignalQuality(device_version=version,scan_ID=run_IDs_1[i],threshold=threshold,channel_name=channel,
                                            local_path=localPath))
idx_allRuns_1 = list(range(0,len(run_IDs_1)))
idx_rejected_1 = score_runs.index('rejected')
idx_accepted_1 = idx_allRuns_1[:idx_rejected_1]
rejected_run_IDs = run_IDs_1[idx_rejected_1]