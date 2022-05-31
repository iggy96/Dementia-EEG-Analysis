# test pipeline for the bruyere dataset

from fn_cfg import *
"""
version = 1.1
filename = '0_1_12072018_1206'
localPath = '/Users/joshuaighalo/Downloads/raw'
y_lim = 1000
figsize = (8,8)
plot_color = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd']
dispIMG = True
titles = ['Fz','Cz','Pz','P07','OZ']

device = importFile.neurocatch()
fileObjects = device.init(version,filename,localPath)
rawEEG = fileObjects[0]
rawEOG = fileObjects[1]
rawEEGEOG = fileObjects[2]
time_period = fileObjects[3]
trigChannel = fileObjects[4]

plots(time_period,rawEEGEOG,titles,figsize,plot_color)
"""



# debugger
#  import raw file 1
localPath = '/Volumes/Backup Plus/EEG_Datasets/laurel_place/dataset'
filename = '0002_2_12122019_1225'
data = filename
localPath = localPath.replace(os.sep, '/')   
localPath = localPath + '/'  
path = localPath+data
os.chdir(path)
for file in (os.listdir(path))[0:2]:
    if file.endswith(".txt"):
        file_path = f"{path}/{file}"
        
with open(file_path, 'r') as f:
    lines = f.readlines(200)

def two_digits(data):
    # electrodes with Nan kOhms 
    reject_impedance = 1000
    fullstring = data
    substring = "NaN"
    if substring in fullstring:
        val3 = reject_impedance
    # electrodes with numeric kOhms
    else:
        val1 = data
        val2 = val1[4:] # delete unnecessary characters from the front
        val2 = val2[:2]
        val3 = int(float(val2)) # delete unnecessary characters from the back then convert to integer
        import math
        # check 1
        digits = int(math.log10(val3))+1 # check number of digits in result
        if digits == 2: # expected result
            val3 = val3
        if digits < 2: # unexpected result 1
            val5 = val1
            val5 = val5[4:]
            val5 = val5[:3]
            val3 = int(float(val5))
        # check 2
        digits = int(math.log10(val3))+1 # check number of digits in result
        if digits == 2: # expected result
            val3 = val3
        if digits < 1: # unexpected result 1
            val6 = val1
            val6 = val6[4:]
            val6 = val6[:4]
            val3 = int(float(val6))
    return val3
    
def three_digits(data):
    # electrodes with Nan kOhms 
    reject_impedance = 1000
    fullstring = data
    substring = "NaN"
    if substring in fullstring:
        val3 = reject_impedance
    # electrodes with numeric kOhms
    else:
        val1 = data
        val2 = val1[4:]
        val2 = val2[:3]
        val3 = int(float(val2))
        import math
        # check 1
        digits = int(math.log10(val3))+1 # check number of digits in result
        if digits == 3: # expected result
            val3 = val3
        if digits < 3: # unexpected result 1
            val5 = val1
            val5 = val5[4:]
            val5 = val5[:4]
            val3 = int(float(val5))
        # check 2
        digits = int(math.log10(val3))+1 # check number of digits in result
        if digits == 3: # expected result
            val3 = val3
        if digits < 2: # unexpected result 1
            val6 = val1
            val6 = val6[4:]
            val6 = val6[:5]
            val3 = int(float(val6))
    return val3


p3 = lines[9] # extract line
print(len(p3))
if len(p3)<=15:
    p3 = two_digits(p3)
else:
    p3 = three_digits(p3)
p4 = lines[11]
if len(p4)<=15:
    p4 = two_digits(p4)
else:
    p4 = three_digits(p4)

p07 = lines[12]
if len(p07)<=15:
    p07 = two_digits(p07)
else:
    p07 = three_digits(p07)

p08 = lines[13]
if len(p08)<=15:
    p08 = two_digits(p08)
else:
    p08 = three_digits(p08)

oz = lines[14]
if len(oz)<=15:
    oz = two_digits(oz)
else:
    oz = three_digits(oz)
    
# extract from metadata file channels that collected eeg data
device_chans = ['FZ','CZ','PZ','PO7','OZ','P3','P4','PO8']
def lcontains(needle_s, haystack_l):
    try: return [i for i in haystack_l if needle_s in i][0]
    except IndexError: return None
metadata_chans = [lcontains(device_chans[i],lines) for i in range(len(device_chans))]

# open meta data file
metaData = open(file_path)
metaData = json.load(metaData)
for i in metaData:
    print(i)
print(metaData['version'])
metaData_chans = metaData['channels']
metaData_imp = metaData['impedances']
# p3,p4,p07,p08,oz
p3 = (metaData_imp[0])['P3']
p4 = (metaData_imp[0])['P4']
p07 = (metaData_imp[0])['PO7']
p08 = (metaData_imp[0])['PO8']
oz = (metaData_imp[0])['OZ']

#  import raw file 1
pathBin = [path+'/'+filename+'.bin']
filenames = glob.glob(pathBin[0])
print(filenames)
data = [np.fromfile(f, dtype=np.float32) for f in filenames]
data_len = len(data)
dataView = data[0]
data1 = data[0]
dataCols = len(metaData_chans)
dataRows = int(len(data1)/dataCols)           
data1 = data1.reshape(dataRows, dataCols)

