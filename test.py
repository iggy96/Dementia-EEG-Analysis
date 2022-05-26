# test pipeline for the bruyere dataset

from fn_cfg import *
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




# debugger
#  import raw file 1
localPath = '/Users/joshuaighalo/Downloads/raw'
filename = '0_1_12072018_1206'
localPath = localPath.replace(os.sep, '/')   
localPath = localPath + '/'  
path = localPath+filename
os.chdir(path)
for file in os.listdir():
    if file.endswith(".json"):
        file_path = f"{path}/{file}"

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