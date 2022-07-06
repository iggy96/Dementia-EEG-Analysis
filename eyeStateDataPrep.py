"""
Generation of eye open/close dataset
Eyes Open and Eyes Close Activity Recognition Using EEG Signals:
For activity recognition, a freely online available EEG-based motor movement and imaginary dataset provided by PhysioNet BCI [12] has been used. 
EEG signals were recorded from all 64 channels for 1 min. 
Two baseline tasks, eyes open (EO) and eyes close (EC) resting state have been used to collect the data from 109 users and are considered in this work. 
In order to detect each activity accurately, the EEG data has been segmented into 10 s. 
Thus, a total of 1308 EEG files (i.e. 654 for EC and 654 for EO) have been created for analysis
link to paper: https://link-springer-com.proxy.lib.sfu.ca/content/pdf/10.1007/978-981-10-9059-2.pdf
link to dataset: https://physionet.org/content/eegmmidb/1.0.0/
"""



from fn_cfg import *
import params as cfg


#   Functions

def extractEDF(file_name,local_directory):
    subfolder = file_name[:-7]
    directory = local_directory + '/' + subfolder + '/' + file_name
    edf_file = mne.io.read_raw_edf(directory)
    raw_data = edf_file.get_data()
    info = edf_file.info
    channelNames = info['ch_names']
    fs = int(info['sfreq'])
    chans_data = raw_data
    Ts = (np.arange(0,len(chans_data.T)/fs,1/fs)).reshape(len(np.arange(0,len(chans_data.T)/fs,1/fs)),1)
    return info,fs,Ts,channelNames,chans_data

def rollingwindow(array,window_size,freq):
    #   Inputs  :   array    - 2D numpy array (d0 = samples, d1 = channels) of filtered EEG data
    #               window_size - size of window to be used for sliding
    #               freq   - step size for sliding window 
    #   Output  :   3D array (columns of array,no of windows,window size)
    def roll_window(array, window_size,freq):
        array = np.array(array,dtype='object')
        shape = (array.shape[0] - window_size + 1, window_size)
        strides = (array.strides[0],) + array.strides
        rolled = np.lib.stride_tricks.as_strided(array, shape=shape, strides=strides)
        return rolled[np.arange(0,shape[0],freq)]
    out_final = []
    for i in range(len(array)):
        out_final.append(roll_window(array[i],window_size,freq))
    out_final = np.asarray(out_final).T
    out_final = out_final.transpose()
    return out_final
    
def featureExtraction(data):
    filtering = filters()
    delta_rhythms = filtering.butterBandPass(data,0.5,4,fs)
    theta_rhythms = filtering.butterBandPass(data,4,8,fs)
    alpha_rhythms = filtering.butterBandPass(data,8,13,fs)
    beta_rhythms = filtering.butterBandPass(data,13,32,fs)
    gamma_rhythms = filtering.butterBandPass(data,32,79,fs)

    std_delta_rhythms = np.std(delta_rhythms,axis=0)
    std_theta_rhythms = np.std(theta_rhythms,axis=0)
    std_alpha_rhythms = np.std(alpha_rhythms,axis=0)
    std_beta_rhythms = np.std(beta_rhythms,axis=0)
    std_gamma_rhythms = np.std(gamma_rhythms,axis=0)

    rms_delta_rhythms = np.sqrt(np.mean(delta_rhythms**2))
    rms_theta_rhythms = np.sqrt(np.mean(theta_rhythms**2))
    rms_alpha_rhythms = np.sqrt(np.mean(alpha_rhythms**2))
    rms_beta_rhythms = np.sqrt(np.mean(beta_rhythms**2))
    rms_gamma_rhythms = np.sqrt(np.mean(gamma_rhythms**2))

    var_delta_rhythms = np.var(delta_rhythms,axis=0)
    var_theta_rhythms = np.var(theta_rhythms,axis=0)
    var_alpha_rhythms = np.var(alpha_rhythms,axis=0)
    var_beta_rhythms = np.var(beta_rhythms,axis=0)
    var_gamma_rhythms = np.var(gamma_rhythms,axis=0)

    kurtosis_delta_rhythms = stats.kurtosis(delta_rhythms,axis=0)
    kurtosis_theta_rhythms = stats.kurtosis(theta_rhythms,axis=0)
    kurtosis_alpha_rhythms = stats.kurtosis(alpha_rhythms,axis=0)
    kurtosis_beta_rhythms = stats.kurtosis(beta_rhythms,axis=0)
    kurtosis_gamma_rhythms = stats.kurtosis(gamma_rhythms,axis=0)
    #features = np.vstack((std_delta_rhythms,std_theta_rhythms,std_alpha_rhythms,std_beta_rhythms,std_gamma_rhythms,rms_delta_rhythms,rms_theta_rhythms,rms_alpha_rhythms,rms_beta_rhythms,rms_gamma_rhythms,var_delta_rhythms,var_theta_rhythms,var_alpha_rhythms,var_beta_rhythms,var_gamma_rhythms,kurtosis_delta_rhythms,kurtosis_theta_rhythms,kurtosis_alpha_rhythms,kurtosis_beta_rhythms,kurtosis_gamma_rhythms)).T
    features = [std_delta_rhythms,std_theta_rhythms,std_alpha_rhythms,std_beta_rhythms,std_gamma_rhythms,rms_delta_rhythms,rms_theta_rhythms,rms_alpha_rhythms,rms_beta_rhythms,rms_gamma_rhythms,var_delta_rhythms,var_theta_rhythms,var_alpha_rhythms,var_beta_rhythms,var_gamma_rhythms,kurtosis_delta_rhythms,kurtosis_theta_rhythms,kurtosis_alpha_rhythms,kurtosis_beta_rhythms,kurtosis_gamma_rhythms]
    return features

def windowEEG(filename,local_dir):
    data = extractEDF(filename,local_dir)
    eeg_data = data[4]
    win_eeg_data = rollingwindow(eeg_data,3200,3200)
    new_eeg = win_eeg_data.reshape(win_eeg_data.shape[0]*win_eeg_data.shape[1],win_eeg_data.shape[2])
    return new_eeg



local_dir = '/Users/joshuaighalo/Downloads/files-2'


files_EO = []
files_EC = []
for root, dirs, files in os.walk(local_dir):
    if files:
        files_EO.append(sorted(files)[0])
        files_EC.append(sorted(files)[2])
files_EO = files_EO[1:]
files_EC = files_EC[1:]

allEEG_EO = []
allEEG_EC = []
for i in range(len(files_EO)):
    allEEG_EO.append(windowEEG(files_EO[i],local_dir))
    allEEG_EC.append(windowEEG(files_EC[i],local_dir))
allEEG_EO = np.array(allEEG_EO)
allEEG_EC = np.array(allEEG_EC)
allEEG_EO = allEEG_EO.reshape(allEEG_EO.shape[0]*allEEG_EO.shape[1],allEEG_EO.shape[2])
allEEG_EC = allEEG_EC.reshape(allEEG_EC.shape[0]*allEEG_EC.shape[1],allEEG_EC.shape[2])



features_EO = []
features_EC = []
for i in range(len(allEEG_EO)):
    features_EO.append(featureExtraction(allEEG_EO[i]))
    features_EC.append(featureExtraction(allEEG_EC[i]))
features_EO = np.array(features_EO)
features_EC = np.array(features_EC)



ones = np.ones(len(features_EO))
ones = ones.reshape(ones.shape[0],1)
zeros = np.zeros(len(features_EC))
zeros = zeros.reshape(zeros.shape[0],1)
features_EO = np.hstack((features_EO,ones))
features_EC = np.hstack((features_EC,zeros))
features = np.vstack((features_EO,features_EC))


df = pd.DataFrame({'std_delta':features[:,0],'std_theta':features[:,1],'std_alpha':features[:,2],'std_beta':features[:,3],'std_gamma':features[:,4],'rms_delta':features[:,5],'rms_theta':features[:,6],'rms_alpha':features[:,7],'rms_beta':features[:,8],'rms_gamma':features[:,9],'var_delta':features[:,10],'var_theta':features[:,11],'var_alpha':features[:,12],'var_beta':features[:,13],'var_gamma':features[:,14],'kurtosis_delta':features[:,15],'kurtosis_theta':features[:,16],'kurtosis_alpha':features[:,17],'kurtosis_beta':features[:,18],'kurtosis_gamma':features[:,19],'eye_status':features[:,20]})
df.to_csv('/Users/joshuaighalo/Documents/BrainNet/Projects/Workspace/results/laurel place/ml_dataset/eye_open_close_windowed_dataset.csv',index=False)