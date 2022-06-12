# test pipeline for the bruyere dataset
from fn_cfg import *
import params as cfg

#localPath = '/Volumes/Backup Plus/EEG_Datasets/laurel_place/dataset'
#filename = '0002_2_12122019_1225'
#version = 1.0

version = 1.1
filename = '0_1_12072018_1206'
localPath = '/Users/joshuaighalo/Downloads/raw'

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
    return np.array(arr,dtype=object).T

def swt(x,wavelet):
    def swt_chans(x):
        coeffs = pywt.swt(x,wavelet,level=10,start_level=1,trim_approx=True)
        return coeffs
    arr = []
    for i in range(len(x.T)):
        arr.append(swt_chans(x[:,i]))
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
    arr_approx = arr_approx
    coefs = arr_detail
    (coefs).insert(0,arr_approx)
    return coefs

def inv_dwt(coeffs,wavelet):
    def inverse_dwt(coeffs,wavelet):
        return waverec(coeffs,wavelet)
    arr = []
    for i in range(len(np.array(coeffs,dtype=object).T)):
        arr.append(inverse_dwt(list(np.array(coeffs,dtype=object)[:,i]),wavelet))
    return  (np.array(arr).T)[:-1,:]


wavelet = 'bior4.4'
coeffs = dwt(notchFilterOutput,wavelet)

threshold_global = global_threshold(notchFilterOutput,coeffs)
threshold_std = std_threshold(coeffs)
coeffs_global = apply_threshold(coeffs,threshold_global)
coeffs_std = apply_threshold(coeffs,threshold_std)
new_signal_global = inv_dwt(coeffs_global,wavelet)
new_signal_std = inv_dwt(coeffs_std,wavelet)

plt.plot(time,new_signal_global[0,:])
plt.show()
plt.plot(time,new_signal_std[0,:])
plt.show()
)

#   Comparative Study of Wavelet-Based Unsupervised Ocular Artifact Removal Techniques for Single-Channel EEG Data
#   Signal to Artifact Ratio (SAR) is a quantification method to measure the amount of artifact removal 
#   in a specific signal after processing with an algorithm [40].
#   SAR is a measure of the amount of artifact removal in a signal.
#   x = EEG signal containing artifact
#   y = EEG signal obtained after running an artifact free algorithm

def sar(x,y):
    return 10*np.log10((np.std(x))/(np.std(x-y)))
def mse(x,y):
    return mean_squared_error(x,y)

sar_global = sar(notchFilterOutput[:,0],new_signal_global[:,0])
sar_std = sar(notchFilterOutput[:,0],new_signal_std[:,0])
print("SAR for global threshold: ",sar_global)
print("SAR for standard deviation threshold: ",sar_std)

mse_global = mse(notchFilterOutput[:,0],new_signal_global[:,0])
mse_std = mse(notchFilterOutput[:,0],new_signal_std[:,0])
print("MSE for global threshold: ",mse_global)
print("MSE for standard deviation threshold: ",mse_std)