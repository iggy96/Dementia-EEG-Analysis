from fn_cfg import *
import params as cfg
from scipy.stats import skew
from sklearn.decomposition import FastICA
import antropy as ant

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

def sample_entropy(data):
    def params(data):
        return ant.sample_entropy(data)
    sample_entropy = []
    for i in range(len(data)):
        sample_entropy.append(params(data[i]))
    return sample_entropy

def confidenceInterval(samples):
#   At 95% significance level, tN -1 = 2.201
    means = np.mean(samples)
    std_dev = np.std(samples)
    standard_error = std_dev/np.sqrt(len(samples))
    lower_95_perc_bound = means - 2.201*standard_error
    upper_95_perc_bound = means + 2.201*standard_error
    return upper_95_perc_bound

def dwt(data,wavelet):
    # A5,D5,D4,D3,D2,D1 = 0-4Hz,8-16Hz,16-32Hz,32-64Hz,64-128Hz
    # A5 is eliminated: OAs are associated with freq ranged < 5Hz
    # D2,D1 are eliminated: MAs are associated with freq ranged > 32Hz
    coeffs = wavedec(data,wavelet,level=5)
    coeffs = np.array(coeffs,dtype=object)
    A5,D5,D4,D3,D2,D1 = coeffs
    zA5,zD2,zD1 = np.zeros((len(A5))),np.zeros(len(D2)),np.zeros((len(D1)))
    value = (np.median(abs(D1))/0.6745)*(np.sqrt(2*np.log(len(data))))
    tD5 = pywt.threshold(D5,value,'soft',substitute=0)
    tD4 = pywt.threshold(D4,value,'soft',substitute=0)
    tD3 = pywt.threshold(D3,value,'soft',substitute=0)
    coeffs = [zA5,tD5,tD4,tD3,zD2,zD1]
    recovered_signal = pywt.waverec(coeffs,wavelet)
    return recovered_signal



hpEEG = icaHighpass(data=rawEEG,cutoff=1,fs=cfg.fs) 

#   Computing ICA components
ica = FastICA(n_components=3, random_state=0, tol=0.0001)
comps = ica.fit_transform(hpEEG)

comps_1_win = rolling_window(comps[:,0],time,1,1)
comps_2_win = rolling_window(comps[:,1],time,1,1)
comps_3_win = rolling_window(comps[:,2],time,1,1)


#   Computing kurtosis of ICA weights
comps_1_kurtosis = kurtosis(comps_1_win,axis=1)
comps_2_kurtosis = kurtosis(comps_2_win,axis=1)
comps_3_kurtosis = kurtosis(comps_3_win,axis=1)

#   Computing skewness of ICA weights
comps_1_skew = skew(comps_1_win,axis=1)
comps_2_skew = skew(comps_2_win,axis=1)
comps_3_skew = skew(comps_3_win,axis=1)

#   Computing peak-to-peak of ICA weights
comps_1_pkpk = peaktopeak(comps_1_win)
comps_2_pkpk = peaktopeak(comps_2_win)
comps_3_pkpk = peaktopeak(comps_3_win)

#   Computing sample entropy of ICA weights
comps_1_sampEN = sample_entropy(comps_1_win)
comps_2_sampEN = sample_entropy(comps_2_win)
comps_3_sampEN = sample_entropy(comps_3_win)


#   Computing CI on kurtosis values
threshold_kurt_1 = confidenceInterval(comps_1_kurtosis)
threshold_kurt_2 = confidenceInterval(comps_2_kurtosis)
threshold_kurt_3 = confidenceInterval(comps_3_kurtosis)
threshold_skew_1 = confidenceInterval(comps_1_skew)
threshold_skew_2 = confidenceInterval(comps_2_skew)
threshold_skew_3 = confidenceInterval(comps_3_skew)
threshold_pkpk_1 = confidenceInterval(comps_1_pkpk)
threshold_pkpk_2 = confidenceInterval(comps_2_pkpk)
threshold_pkpk_3 = confidenceInterval(comps_3_pkpk)
threshold_sampEN_1 = confidenceInterval(comps_1_sampEN)
threshold_sampEN_2 = confidenceInterval(comps_2_sampEN)
threshold_sampEN_3 = confidenceInterval(comps_3_sampEN)



#   Extract epochs
bool_ArtfCompsKurt_1 = [comps_1_kurtosis>threshold_kurt_1]
idx_ArtfCompsKurt_1 = np.asarray(np.where(bool_ArtfCompsKurt_1[0]==True))
bool_ArtfCompsKurt_2 = [comps_2_kurtosis>threshold_kurt_2]
idx_ArtfCompsKurt_2 = np.asarray(np.where(bool_ArtfCompsKurt_2[0]==True))
bool_ArtfCompsKurt_3 = [comps_3_kurtosis>threshold_kurt_3]
idx_ArtfCompsKurt_3 = np.asarray(np.where(bool_ArtfCompsKurt_3[0]==True))

bool_ArtfCompsSkew_1 = [comps_1_skew>threshold_skew_1]
idx_ArtfCompsSkew_1 = np.asarray(np.where(bool_ArtfCompsSkew_1[0]==True))
bool_ArtfCompsSkew_2 = [comps_2_skew>threshold_skew_2]
idx_ArtfCompsSkew_2 = np.asarray(np.where(bool_ArtfCompsSkew_2[0]==True))
bool_ArtfCompsSkew_3 = [comps_3_skew>threshold_skew_3]
idx_ArtfCompsSkew_3 = np.asarray(np.where(bool_ArtfCompsSkew_3[0]==True))

bool_ArtfCompsPkPk_1 = [comps_1_pkpk>threshold_pkpk_1]
idx_ArtfCompsPkPk_1 = np.asarray(np.where(bool_ArtfCompsPkPk_1[0]==True))
bool_ArtfCompsPkPk_2 = [comps_2_pkpk>threshold_pkpk_2]
idx_ArtfCompsPkPk_2 = np.asarray(np.where(bool_ArtfCompsPkPk_2[0]==True))
bool_ArtfCompsPkPk_3 = [comps_3_pkpk>threshold_pkpk_3]
idx_ArtfCompsPkPk_3 = np.asarray(np.where(bool_ArtfCompsPkPk_3[0]==True))

bool_ArtfCompsSampEN_1 = [comps_1_sampEN>threshold_sampEN_1]
idx_ArtfCompsSampEN_1 = np.asarray(np.where(bool_ArtfCompsSampEN_1[0]==True))
bool_ArtfCompsSampEN_2 = [comps_2_sampEN>threshold_sampEN_2]
idx_ArtfCompsSampEN_2 = np.asarray(np.where(bool_ArtfCompsSampEN_2[0]==True))
bool_ArtfCompsSampEN_3 = [comps_3_sampEN>threshold_sampEN_3]
idx_ArtfCompsSampEN_3 = np.asarray(np.where(bool_ArtfCompsSampEN_3[0]==True))



#   merge indices arrays
idx_ArtfComps_1 = np.concatenate((idx_ArtfCompsKurt_1,idx_ArtfCompsSkew_1,idx_ArtfCompsPkPk_1,idx_ArtfCompsSampEN_1),axis=1).T
idx_ArtfComps_1 = np.unique(idx_ArtfComps_1,axis=0)
idx_ArtfComps_2 = np.concatenate((idx_ArtfCompsKurt_2,idx_ArtfCompsSkew_2,idx_ArtfCompsPkPk_2,idx_ArtfCompsSampEN_2),axis=1).T
idx_ArtfComps_2 = np.unique(idx_ArtfComps_2,axis=0)
idx_ArtfComps_3 = np.concatenate((idx_ArtfCompsKurt_3,idx_ArtfCompsSkew_3,idx_ArtfCompsPkPk_3,idx_ArtfCompsSampEN_3),axis=1).T
idx_ArtfComps_3 = np.unique(idx_ArtfComps_3,axis=0)


""" 
Apply Discrete Wavelet Transform to epochs of ICs
remove noisy componenets within epoch
apply inverse dwt to retrieve denoised epochs
"""

wavelet = 'sym9'    
inv_comps_1 = []
inv_comps_2 = []
inv_comps_3 = []

artf_comps_1 = comps_1_win[idx_ArtfComps_1]
artf_comps_1 = artf_comps_1.reshape(artf_comps_1.shape[0],artf_comps_1.shape[1]*artf_comps_1.shape[2])
artf_comps_2 = comps_2_win[idx_ArtfComps_2]
artf_comps_2 = artf_comps_2.reshape(artf_comps_2.shape[0],artf_comps_2.shape[1]*artf_comps_2.shape[2])
artf_comps_3 = comps_3_win[idx_ArtfComps_3]
artf_comps_3 = artf_comps_3.reshape(artf_comps_3.shape[0],artf_comps_3.shape[1]*artf_comps_3.shape[2])

for i in range(len(idx_ArtfComps_1)):
    inv_comps_1.append(dwt(artf_comps_1[i],wavelet))
inv_comps_1 = np.array(inv_comps_1)
for i in range(len(idx_ArtfComps_2)):
    inv_comps_2.append(dwt(artf_comps_2[i],wavelet))
inv_comps_2 = np.array(inv_comps_2)
for i in range(len(idx_ArtfComps_3)):
    inv_comps_3.append(dwt(artf_comps_3[i],wavelet))
inv_comps_3 = np.array(inv_comps_3)

a1 = np.arange(10)
a1  =a1.reshape(2,5).T
a2 = np.arange(20)
a2 = a2[9:19]
a2 = a2.reshape(2,5).T
idx = np.array([0,3,4])

for i in range(len(idx)):
    a2[idx[i],:] = a1[idx[i],:]


for i in range(len(idx_ArtfComps_1)):
    inv_comps_1[i,:] = artf_comps_1[i,:]
