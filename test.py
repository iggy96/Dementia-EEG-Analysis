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


wavelet = 'coif5'
coeffs = dwt(notchFilterOutput,wavelet)
threshold_global = global_threshold(notchFilterOutput,coeffs)
threshold_std = std_threshold(coeffs)
coeffs_global = apply_threshold(coeffs,threshold_global)
coeffs_std = apply_threshold(coeffs,threshold_std)
new_signal_global = inv_dwt(coeffs_global,wavelet)
new_signal_std = inv_dwt(coeffs_std,wavelet)

plt.plot(time,new_signal_global[:,0],color='r',label='Global Threshold')
plt.legend()
plt.show()
plt.plot(time,new_signal_std[:,0],color='b',label='STD Threshold')
plt.legend()
plt.show()

adaptiveFilterOutput = filtering.adaptive(new_signal_global,rawEOG)
plots(time,adaptiveFilterOutput,titles=cfg.channelNames,figsize=cfg.figure_size,pltclr=cfg.plot_color)
bandPassFilterOutput = filtering.butterBandPass(adaptiveFilterOutput,lowcut=cfg.highPass,highcut=cfg.lowPass,fs=cfg.fs)
plots(time,bandPassFilterOutput,titles=cfg.channelNames,figsize=cfg.figure_size,pltclr=cfg.plot_color)

erps = erpExtraction()
N1P3 = erps.N100P300(trigOutput,bandPassFilterOutput,time,stimTrig=cfg.stimTrig,clip=cfg.clip)
N4 = erps.N400(trigOutput,bandPassFilterOutput,time,stimTrig=cfg.stimTrig,clip=cfg.clip)
N1P3_Fz = N1P3[0]
N1P3_Cz = N1P3[1]
N1P3_Pz = N1P3[2]
N4_Fz = N4[0]
N4_Cz = N4[1]
N4_Pz = N4[2]
erp_latency = np.array(np.linspace(start=-100, stop=900, num=len(N1P3_Fz[0]),dtype=object),dtype=object)
std_N1P3_Fz = N1P3_Fz[0]

plot_ERPs(N1P3_Fz[0],N1P3_Fz[1],erp_latency,'N1P3_Fz','Latency (ms)','Amplitude (uV)','std','dev','b','r',10)
plot_ERPs(N1P3_Cz[0],N1P3_Cz[1],erp_latency,'N1P3_Cz','Latency (ms)','Amplitude (uV)','std','dev','b','r',10)
plot_ERPs(N1P3_Pz[0],N1P3_Pz[1],erp_latency,'N1P3_Pz','Latency (ms)','Amplitude (uV)','std','dev','b','r',10)
plot_ERPs(N4_Fz[0],N4_Fz[1],erp_latency,'N4_Fz','Latency (ms)','Amplitude (uV)','std','dev','b','r',10)
plot_ERPs(N4_Cz[0],N4_Cz[1],erp_latency,'N4_Cz','Latency (ms)','Amplitude (uV)','std','dev','b','r',10)
plot_ERPs(N4_Pz[0],N4_Pz[1],erp_latency,'N4_Pz','Latency (ms)','Amplitude (uV)','std','dev','b','r',10)

"""
#   SWT Wavelet Denoising
input_data = notchFilterOutput[:,0]

def length_check_1(input_data):
    if len(input_data)%2 == 1:
        input_data = np.pad(input_data, (0, 1), 'constant')
    else:
        input_data = input_data
    return input_data

input_data = length_check_1(input_data)
coeffs_swt = pywt.swt(input_data,wavelet,level=1,trim_approx=True)
threshold_global_approx = (np.median(abs(coeffs_swt[0]))/0.6745)*(np.sqrt(2*np.log(len(input_data))))
threshold_global_detail = (np.median(abs(coeffs_swt[1]))/0.6745)*(np.sqrt(2*np.log(len(input_data))))
threshold_std_approx = 1.5*np.std(coeffs_swt[0])
threshold_std_detail = 1.5*np.std(coeffs_swt[1])

def swt_ApplyThreshold(coeffs,threshold_approx,threshold_detail):
    def approx_ApplyThreshold(coeffs,threshold_approx):
        coeff_1 = coeffs[0]
        coeff_1[abs(coeff_1)>threshold_approx] = 0
        return list(coeff_1)
    def detail_ApplyThreshold(coeffs,threshold_detail):
        coeff_1 = coeffs[1]
        coeff_1[abs(coeff_1)>threshold_detail] = 0
        return list(coeff_1)
    approx = approx_ApplyThreshold(coeffs,threshold_approx)
    detail = detail_ApplyThreshold(coeffs,threshold_detail)
    coefs = [approx,detail]
    return coefs

coeffs_global = swt_ApplyThreshold(coeffs_swt,threshold_global_approx,threshold_global_detail)
coeffs_std = swt_ApplyThreshold(coeffs_swt,threshold_std_approx,threshold_std_detail)
new_signal_global = pywt.iswt(np.array(coeffs_global),wavelet)
new_signal_global = new_signal_global.reshape(len(new_signal_global),1)
new_signal_std = pywt.iswt(np.array(coeffs_std),wavelet)
new_signal_std = new_signal_std.reshape(len(new_signal_std),1)

def length_check_2(signal):
    if len(signal) % 2 == 0:
        signal = np.delete(signal,-1)
    else:
        pass
    return signal

new_signal_global = length_check_2(new_signal_global)
new_signal_std = length_check_2(new_signal_std)

plt.plot(time,new_signal_global,color='r',label='Global Threshold')
plt.legend()
plt.show()
plt.plot(time,new_signal_std,color='r',label='STD Threshold')
plt.legend()
plt.show()
"""


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