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

def dwt(x,wavelet):
    def dwt_chans(x):
        coeffs = wavedec(x,wavelet,level=10)
        return coeffs
    arr = []
    for i in range(len(x.T)):
        arr.append(dwt_chans(x[:,i]))
    return np.array(arr).T

def global_threshold(data,coeffs):
    def coeffs_approx(data,coeffs):
        return (np.median(abs(coeffs[0]))/0.6745)*(np.sqrt(2*np.log(len(data))))
    def coeffs_detail(data,coeffs):
        return (np.median(abs(coeffs[1]))/0.6745)*(np.sqrt(2*np.log(len(data))))
    arr_approx = [ ]
    arr_detail = [ ]
    for i in range(len(coeffs.T)):
        arr_approx.append(coeffs_approx(data,coeffs[:,i]))
        arr_detail.append(coeffs_detail(data,coeffs[:,i]))
    return np.vstack((arr_approx,arr_detail))

def std_threshold(coeffs):
    def std_approx(coeffs):
        return 1.5*np.std(coeffs[0])
    def std_detail(coeffs):
        return 1.5*np.std(coeffs[1])
    arr_approx = [ ]
    arr_detail = [ ]
    for i in range(len(coeffs.T)):
        arr_approx.append(std_approx(coeffs[:,i]))
        arr_detail.append(std_detail(coeffs[:,i]))
    return np.vstack((arr_approx,arr_detail))

def apply_threshold(coeffs,threshold):
    def apply_threshold_approx(coeffs,threshold):
        coeffs[0][abs(coeffs[0])>threshold[0]] = 0
        coeffs_approx = coeffs[0]
        return coeffs_approx
    def apply_threshold_detail(coeffs,threshold):
        coeffs = coeffs[1:len(coeffs)]
        coeffs[0][abs(coeffs[0])>threshold[1]] = 0
        coeffs[1][abs(coeffs[1])>threshold[1]] = 0
        coeffs[2][abs(coeffs[2])>threshold[1]] = 0
        coeffs[3][abs(coeffs[3])>threshold[1]] = 0
        coeffs[4][abs(coeffs[4])>threshold[1]] = 0
        coeffs[5][abs(coeffs[5])>threshold[1]] = 0
        coeffs[6][abs(coeffs[6])>threshold[1]] = 0
        coeffs[7][abs(coeffs[7])>threshold[1]] = 0
        coeffs[8][abs(coeffs[8])>threshold[1]] = 0
        coeffs[9][abs(coeffs[9])>threshold[1]] = 0
        return coeffs
    arr_approx = [ ]
    arr_detail = [ ]
    for i in range(len(coeffs.T)):
        arr_approx.append(apply_threshold_approx(coeffs[:,i],threshold[:,i]))
    for i in range((len(coeffs.T))):
        arr_detail.append(apply_threshold_detail(coeffs[:,i],threshold[:,i]))
    arr_detail = list(np.array(arr_detail).T)
    arr_approx = list(zip(arr_approx))
    return arr_approx,




wavelet = 'haar'
coeffs = dwt(notchFilterOutput,wavelet)

threshold_global = global_threshold(notchFilterOutput,coeffs)
threshold_std = std_threshold(coeffs)
coeffs_approx_global = (apply_threshold(coeffs,threshold_global))[0]
coeffs_detail_global = (apply_threshold(coeffs,threshold_global))[1]
coeffs_approx_std = (apply_threshold(coeffs,threshold_std))[0]
coeffs_detail_std = (apply_threshold(coeffs,threshold_std))[1]

ctest = list(zip(coeffs_approx_global))
list1 = ctest
list2 = coeffs_detail_global
(list2).insert(0,list1)


"""
def dwt(x,wavelet):
    coeffs = wavedec(x,wavelet,level=10)
    return coeffs

def global_threshold(data,coeffs):
    approx = (np.median(abs(coeffs[0]))/0.6745)*(np.sqrt(2*np.log(len(data))))
    detail = (np.median(abs(coeffs[1]))/0.6745)*(np.sqrt(2*np.log(len(data))))
    return np.vstack((approx,detail))

def std_threshold(coeffs):
    approx = 1.5*np.std(coeffs[0])
    detail =  1.5*np.std(coeffs[1])
    return np.vstack((approx,detail))

def apply_threshold(coeffs,threshold):
    def apply_threshold_approx(coeffs,threshold):
        coeffs[0][abs(coeffs[0])>threshold[0]] = 0
        coeffs_approx = coeffs[0]
        return coeffs_approx

    def apply_threshold_detail(coeffs,threshold):
        coeffs = coeffs[1:len(coeffs)]
        coeffs[abs(coeffs)>threshold[1]] = 0
        return coeffs 

    coeffs_detail = []
    for i in range(len(coeffs)-1):
        coeffs_detail.append(apply_threshold_detail(coeffs[i],threshold))
    coeffs_detail = np.array(coeffs_detail)
    coeffs_approx = apply_threshold_approx(coeffs,threshold)
    return coeffs_approx,coeffs_detail

def inverse_dwt(coeffs,wavelet):
    return waverec(coeffs,wavelet)


wavelet = 'haar'
data = notchFilterOutput[:,0]

coeffs = dwt(data,wavelet)
threshold_global = global_threshold(data,coeffs)
threshold_std = std_threshold(coeffs)
coeffs_approx_global = (apply_threshold(coeffs,threshold_global))[0]
coeffs_detail_global = (apply_threshold(coeffs,threshold_global))[1]
coeffs_approx_std = (apply_threshold(coeffs,threshold_std))[0]
coeffs_detail_std = (apply_threshold(coeffs,threshold_std))[1]


"""

def swt(x):
    def swt_epochs(x):
        coeffs = pywt.swt(x, 'haar', level=2, trim_approx=True)
        return coeffs
    arr = []
    for i in range(len(x)):
        arr.append(swt_epochs(x[i]))
    return arr

epochs_swt = swt(epochs)






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


















