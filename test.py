from fn_cfg import *
import params as cfg
from scipy.stats import skew
from sklearn.decomposition import FastICA

"Laurel Place Dataset"
localPath = '/Users/joshuaighalo/Downloads/brainNet_datasets/laurel_place/cleaned_dataset'
img_save_path = '/Users/joshuaighalo/Documents/BrainNet/Projects/Workspace/results/laurel place/images/single_averages/2.0'
filename = '0002_1_12122019_1219'
version = 1.0
image_name = "Dementia Scan"
cutoff = 10
device = importFile.neurocatch()
fileObjects = device.init(version,filename,localPath,dispIMG=False)
rawEEG = fileObjects[0]
rawEOG = fileObjects[1]
rawEEGEOG = fileObjects[2]
time = fileObjects[3]



#   Implement high pass filter @ 1Hz
def icaHighpass(data,cutoff,fs):
    def params_fnc(data,cutoff,fs,order=4):
        nyq = 0.5 * fs
        normal_cutoff = cutoff / nyq
        b, a = signal.butter(order, normal_cutoff, btype='high', analog=False)
        y = signal.filtfilt(b, a, data)
        return y
    filterEEG = []
    for i in range(len(data.T)):
        filterEEG.append(params_fnc(data.T[i],cutoff,fs))
    filterEEG = np.array(filterEEG).T
    return filterEEG

hpEEG = icaHighpass(data=rawEEG,cutoff=1,fs=cfg.fs) 

#   Computing ICA components
ica = FastICA(n_components=3, random_state=0, tol=0.0001)
comps = ica.fit_transform(hpEEG)

comps_1_win = rolling_window(comps[:,0],time,10,10)
comps_2_win = rolling_window(comps[:,1],time,10,10)
comps_3_win = rolling_window(comps[:,2],time,10,10)


#   Computing kurtosis of ICA weights
comps_1_kurtosis = kurtosis(comps_1_win,axis=1)
comps_2_kurtosis = kurtosis(comps_2_win,axis=1)
comps_3_kurtosis = kurtosis(comps_3_win,axis=1)
comps_1_skew = skew(comps_1_win,axis=1)
comps_2_skew = skew(comps_2_win,axis=1)
comps_3_skew = skew(comps_3_win,axis=1)


#   Computing CI on kurtosis values
def confidenceInterval(samples):
#   At 95% significance level, tN -1 = 2.201
    means = np.mean(samples)
    std_dev = np.std(samples)
    standard_error = std_dev/np.sqrt(len(samples))
    lower_95_perc_bound = means - 2.201*standard_error
    upper_95_perc_bound = means + 2.201*standard_error
    return upper_95_perc_bound


threshold_kurt_1 = confidenceInterval(comps_1_kurtosis)
threshold_kurt_2 = confidenceInterval(comps_2_kurtosis)
threshold_kurt_3 = confidenceInterval(comps_3_kurtosis)
threshold_skew_1 = confidenceInterval(comps_1_skew)
threshold_skew_2 = confidenceInterval(comps_2_skew)
threshold_skew_3 = confidenceInterval(comps_3_skew)

#   