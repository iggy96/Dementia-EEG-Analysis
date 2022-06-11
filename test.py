# test pipeline for the bruyere dataset

from email.mime import base
from bleach import clean
from fn_cfg import *
import params as cfg

#localPath = '/Volumes/Backup Plus/EEG_Datasets/laurel_place/dataset'
#filename = '0002_2_12122019_1225'
#version = 1.0

version = 1.1
filename = '0_1_12072018_1206'
localPath = '/Users/joshuaighalo/Downloads/raw'


dispIMG = True


device = importFile.neurocatch()
fileObjects = device.init(version,filename,localPath)
rawEEG = fileObjects[0]
rawEOG = fileObjects[1]
rawEEGEOG = fileObjects[2]
time = fileObjects[3]
trigOutput = fileObjects[4]
plots(time,rawEEG,titles=cfg.channelNames,figsize=cfg.figure_size,pltclr=cfg.plot_color)
spectogramPlot(rawEEG,fs,nfft=cfg.nfft,nOverlap=cfg.noverlap,figsize=(16,6),subTitles=cfg.channelNames,title='Music Therapy Group 11')
filtering = filters()
notchFilterOutput = filtering.notch(rawEEG,line,fs)
plots(time,notchFilterOutput,titles=cfg.channelNames,figsize=cfg.figure_size,pltclr=cfg.plot_color)
spectogramPlot(notchFilterOutput,fs,nfft=cfg.nfft,nOverlap=cfg.noverlap,figsize=(16,6),subTitles=cfg.channelNames,title='Music Therapy Group 11')









"""
Probability Mapping Based Artifact Detection and Wavelet Denoising based 
Artifact Removal from Scalp EEG for BCI Applications
"""
filtering = filters()
notchFilterOutput = filtering.notch(rawEEG,line,fs)
epochs = rolling_window(notchFilterOutput[:,0],time,1,1)

#   Artifact detection

def entropy(data):
    arr = []
    for i in range(len(data)):
        arr.append(float(stats.entropy(data[i],base=2)))
    return np.array(arr)

def kurtosis(data):
    arr = []
    for i in range(len(data)):
        arr.append(scipy.stats.kurtosis(data[i]))
    return np.array(arr)

def skewness(data):
    arr = []
    for i in range(len(data)):
        arr.append(scipy.stats.skew(data[i]))
    return np.array(arr)

epochs_entropy = entropy(epochs)
epochs_kurtosis = kurtosis(epochs)
epochs_skewness = skewness(epochs)

time_t = np.arange(0, 1, 1/500)
test_clean = 100 * np.sin(2 * np.pi * 2 * time_t + 0)
test_noisy = 500 * np.sin(2 * np.pi * 2 * time_t + 0)
plt.plot(time_t, test_clean, 'b-', label='clean')
plt.plot(time_t, test_noisy, 'r-', label='noisy')
plt.legend()
plt.show()

print("entropy clean: ", entropy(test_clean.reshape(1,len(test_clean))))
print("entropy noisy: ", entropy(test_noisy.reshape(1,len(test_noisy))))










"""
# extract clean eeg from original signal
def length(x):
    arr1 = []
    for i in range(len(x)):
        arr1.append(len(x[i]))
    return arr1
def index(x,y):
    arr1 = []
    for i in range(len(y)):
        arr1.append(x[y[i]])
    return np.array(arr1)
def median(data):
    arr = []
    for i in range(len(data)):
        arr.append(np.median(data[i]))
    return np.array(arr)
def indexNearestCleanEEG(source,bank):
    def nearestIndex(num,array):
        return np.where(abs(array-num)==abs(array-num).min())[0]
    arr = []
    for i in range(len(source)):
        arr.append(nearestIndex(source[i],bank))
    return np.array(arr)


origEEG = rawEEG[:,0]
origEEG[abs(origEEG) <= 50] = np.nan
idx_cleanEEG = np.argwhere(np.isnan(origEEG))
idx_groupCleanEEG = np.split(idx_cleanEEG, np.cumsum( np.where(idx_cleanEEG[1:] - idx_cleanEEG[:-1] > 1) )+1)
idx_NonZeroGroup = np.where(np.array(length(idx_groupCleanEEG)) > 0)[0].ravel()
idx_groupCleanEEG = index(idx_groupCleanEEG,idx_NonZeroGroup)
len_group = length(idx_groupCleanEEG)
len_max = np.amax(len_group)
idx_len_max = np.argwhere(len_group == len_max)
ext_idx = idx_groupCleanEEG[idx_len_max[0][0]]

median_idx_groupCleanEEG = median(idx_groupCleanEEG)

device = importFile.neurocatch()
fileObjects = device.init(version,filename,localPath)
rawEEG = fileObjects[0]
segment_CleanEEG = rawEEG[:,0][ext_idx]

plt.plot(segment_CleanEEG)
plt.show()


#   extract artifacted eeg from original signal
origEEG = rawEEG[:,0]
origEEG[abs(origEEG) > 50] = np.nan
idx_ArtEEG = np.argwhere(np.isnan(origEEG))
idx_groupArtEEG = np.split(idx_ArtEEG, np.cumsum( np.where(idx_ArtEEG[1:] - idx_ArtEEG[:-1] > 1) )+1)
len_group = np.array(length(idx_groupArtEEG))
idx_NonZeroGroup = np.where(len_group > 0)[0].ravel()
idx_groupArtEEG = index(idx_groupArtEEG,idx_NonZeroGroup)

median_idx_groupArtEEG = median(idx_groupArtEEG)

device = importFile.neurocatch()
fileObjects = device.init(version,filename,localPath)
rawEEG = fileObjects[0]
segArtEEG = index(rawEEG[:,0],idx_groupArtEEG)

plt.plot(segArtEEG[0])
plt.show()

#   evaluate the reference clean signals to be used
#   use the indices of clean eeg close to the artifact indices 
#   to extract the clean eeg signals
idx_nearest_cleanEEG = indexNearestCleanEEG(median_idx_groupArtEEG,median_idx_groupCleanEEG)
idx_nearest_cleanEEG = index(idx_groupCleanEEG,idx_nearest_cleanEEG)

device = importFile.neurocatch()
fileObjects = device.init(version,filename,localPath)
rawEEG = fileObjects[0]

refCleanEEG =  index(rawEEG[:,0],idx_nearest_cleanEEG.ravel())



#   wavelets decompositions

original_EEG = rawEEG[:,2]
clean_cA,clean_cD5,clean_cD4,clean_cD3,clean_cD2,clean_cD1 = wavedec(refCleanEEG[0],'bior4.4',level=5)
art_cA,art_cD5,art_cD4,art_cD3,art_cD2,art_cD1 = wavedec(segArtEEG[0],'bior4.4',level=5)


cdf_clean_cD5 = np.cumsum(clean_cD5)
cdf_clean_cD4 = np.cumsum(clean_cD4)
cdf_clean_cD3 = np.cumsum(clean_cD3)
cdf_clean_cD2 = np.cumsum(clean_cD2)
cdf_clean_cD1 = np.cumsum(clean_cD1)
cdf_art_cD5 = np.cumsum(art_cD5)
cdf_art_cD4 = np.cumsum(art_cD4)
cdf_art_cD3 = np.cumsum(art_cD3)
cdf_art_cD2 = np.cumsum(art_cD2)
cdf_art_cD1 = np.cumsum(art_cD1)

rec_cdf_clean_cD5 = np.reciprocal(cdf_clean_cD5)
rec_cdf_clean_cD4 = np.reciprocal(cdf_clean_cD4)
rec_cdf_clean_cD3 = np.reciprocal(cdf_clean_cD3)
rec_cdf_clean_cD2 = np.reciprocal(cdf_clean_cD2)
rec_cdf_clean_cD1 = np.reciprocal(cdf_clean_cD1)

tm_cD5 = np.multiply(rec_cdf_clean_cD5,cdf_art_cD5)
tm_cD4 = np.multiply(rec_cdf_clean_cD4,cdf_art_cD4)
tm_cD3 = np.multiply(rec_cdf_clean_cD3,cdf_art_cD3)
tm_cD2 = np.multiply(rec_cdf_clean_cD2,cdf_art_cD2)
tm_cD1 = np.multiply(rec_cdf_clean_cD1,cdf_art_cD1)
"""


















