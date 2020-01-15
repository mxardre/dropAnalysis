
# coding: utf-8

# In[ ]:




# In[2]:

#import the droplet file and load them in a dataframe
import pandas as pd
import numpy as np
from os import listdir
from os.path import isfile, join
import matplotlib.pyplot as plt

dt = np.dtype([('time','<i4'),]+[(x, '<i2') for x in ['pmt1','pmt2','pmt3']])

#path='/Users/maxime/Documents/experiences/milli/MILLIDROP_manip/2017-10-11_CAABipy/analysis/'
#path='/Users/maxime/Documents/experiences/milli/MILLIDROP_manip/2017-10-13_CAABipy_WT/analysis/'
#path='/Users/maxime/Documents/experiences/milli/MILLIDROP_manip/2017-10-17_CAABipy_pvdS/analysis/'
#path='/Users/maxime/Documents/experiences/milli/MILLIDROP_manip/2017-10-17_CAABipy_pvdS/analysis/'
path='/Users/maxime/Documents/experiences/milli/MILLIDROP_manip/2017-10-19_dilution_CAA_WT/analysis/'

drpfiles=[path+'droplets/' + f for f in listdir(path + 'droplets/') if isfile(join(path + 'droplets/',f))]

i=0
df={}
for file in drpfiles:
    df[i]=pd.read_csv(file)
    i+=1

    #create the template file in hiccup mode
tpfile=pd.read_csv(path+'template.csv', header=4)
tpFileOrd=tpfile.set_index('order')

dropMap=[]

for i in range(0,len(tpFileOrd)-1,2):
    for j in range(2*(tpFileOrd.droplet_number[i])):
        if j%2==0 :
            dropMap.append([tpFileOrd.well[i], tpFileOrd.description[i]])
        else : 
            dropMap.append([tpFileOrd.well[i+1], tpFileOrd.description[i+1]])

dropMap=np.array(dropMap)
#dropMap[:,1]
#dropMap = np.delete(dropMap, (0), axis=0)
#print(dropMap[0:25])


# In[ ]:




# In[3]:

#plot all the data roughly
import matplotlib.pyplot as plt

for i in range(0,len(df)):
    p=df[i]
    plt.scatter(p.time/3600, p.fluo_3_mean)

plt.show()


# In[ ]:




# In[6]:

def plotLabel(df, dropMap, label, channel,low, high):

    for i,content in enumerate(dropMap[:,1]):
        if content==label:
            if channel==2:
                p=df[i]
                plt.scatter(p.time/3600, p.fluo_2_mean)
            if channel==3:
                p=df[i]
                plt.scatter(p.time/3600, p.fluo_3_mean)
            if channel==1:
                p=df[i]
                plt.scatter(p.time/3600, p.fluo_1_mean)
    plt.ylim([low, high])
    #plt.xlim([0, 48])
    plt.title(label+' channel='+str(channel))
    plt.show()


# In[8]:

#label='PvdS-25'
#label='PvdS-5'
#label='Empty'
#label='KB'
#label='SBW25-25'

#EXP20171004_1537
#CAA1mM
#Empty 
#CAA100uM 
#CAA0M 
#M9-glu1mM 
#M9-glu100uM 
#M9-glu0M 
#pvd
label =['CAA1mM', 'CAA100uM', 'CAA0M', 'Empty', ',CAA400uM', ',CAA700uM', 'CAA0M-eth', 'CAA0M']
channel=3
for tmp in label:
    plotLabel(df, dropMap, tmp,channel, 0, 5)

channel=1
for tmp in label:
    plotLabel(df, dropMap, tmp,channel, 0, 5)
    
channel=2
for tmp in label:
    plotLabel(df, dropMap, tmp,channel, 0, 5)


# In[9]:

import matplotlib.pyplot as plt
def plotLabelSeveral(df, dropMap, label, channel,low=0, high=5, namePath=''):
    choice=['blue','green','red','cyan','magenta','yellow','black']

    colors=choice[0:len(label)];
    
    for i,content in enumerate(dropMap[:,1]):
        for j,tmp in enumerate(label):
            if content==tmp:
                if channel=='mCherry':
                    p=df[i]
                    plt.scatter(p.time/3600, p.fluo_2_mean,c=colors[j], cmap=plt.cm.RdYlGn)
                if channel=='GFP':
                    p=df[i]
                    plt.scatter(p.time/3600, p.fluo_3_mean,c=colors[j], cmap=plt.cm.RdYlGn)
                if channel=='PVD':
                    p=df[i]
                    plt.scatter(p.time/3600, p.fluo_1_mean,c=colors[j], cmap=plt.cm.RdYlGn)

    fig1 = plt.gcf()
    plt.ylim([low, high])
    #plt.xlim([0, 24])
    plt.title(' '.join(label)+' '+ ' '.join(choice[0:len(label)])+"\nchannel="+ str(channel) )
    plt.xlabel('time (h)')
    plt.ylabel('fluo (ua)')
    plt.show()
    plt.draw()
    if namePath!='':
        fig1.savefig(namePath)


# In[10]:

#label =['CAA1mM', 'CAA100uM', 'CAA0M', 'Empty', 'M9-glu1mM', 'M9-glu100uM', 'M9-glu0M', 'pvd']
label =['CAA1mM', 'CAA100uM', 'CAA0M', 'Empty', ',CAA400uM', ',CAA700uM', 'CAA0M-eth', 'CAA0M']

channel='mCherry'
plotLabelSeveral(df, dropMap, ['CAA1mM', 'CAA100uM','CAA0M'],channel, 0, 5,path+'CAAChannelmC.jpg')
#plotLabelSeveral(df, dropMap, ['M9-glu1mM', 'M9-glu100uM', 'M9-glu0M'],channel, 0, 5,path+'M9gluChannelmC.jpg')
plotLabelSeveral(df, dropMap, ['pvd'],channel, 0, 5,path+'pvdChannelmC.jpg')

channel='PVD'
plotLabelSeveral(df, dropMap, ['CAA1mM', 'CAA100uM','CAA0M'],channel, 0, 5,path+'CAAChannelPVD.jpg')
#plotLabelSeveral(df, dropMap, ['M9-glu1mM', 'M9-glu100uM', 'M9-glu0M'],channel, 0, 5,path+'M9gluChannelPVD.jpg')
plotLabelSeveral(df, dropMap, ['pvd'],channel, 0, 5,path+'pvdChannelPVD.jpg')

channel='GFP'
plotLabelSeveral(df, dropMap, ['CAA1mM', 'CAA100uM','CAA0M'],channel, 0, 5,path+'CAAChannelGFP.jpg')
#plotLabelSeveral(df, dropMap, ['M9-glu1mM', 'M9-glu100uM', 'M9-glu0M'],channel, 0, 5,path+'M9gluChannelGFP.jpg')
plotLabelSeveral(df, dropMap, ['pvd'],channel, 0, 5,path+'pvdChannelGFP.jpg')


# In[54]:

import matplotlib.pyplot as plt
def plotLabelSeveralLog(df, dropMap, label, channel,low=0, high=5, namePath=''):
    choice=['blue','green','red','cyan','magenta','yellow','black']

    colors=choice[0:len(label)];
    
    for i,content in enumerate(dropMap[:,1]):
        for j,tmp in enumerate(label):
            if content==tmp:
                if channel=='mCherry':
                    p=df[i]
                    plt.scatter(p.time/3600[::2], np.log(p.fluo_2_mean[::2]),c=colors[j], cmap=plt.cm.RdYlGn)
                if channel=='GFP':
                    p=df[i]
                    plt.scatter(p.time[::2]/3600, np.log(p.fluo_3_mean[::2]),c='blue', cmap=plt.cm.RdYlGn)
                    #plt.scatter(p.time/3600, np.log(p.fluo_3_mean),c='blue', cmap=plt.cm.RdYlGn)
                    plt.plot(p.time/3600, np.log(p.fluo_3_mean))#,c='red', cmap=plt.cm.RdYlGn)
   
                if channel=='PVD':
                    p=df[i]
                    plt.scatter(p.time/3600, np.log(p.fluo_1_mean),c=colors[j], cmap=plt.cm.RdYlGn)

    fig1 = plt.gcf()
    plt.ylim([np.log(low), np.log(high)])
    #plt.xlim([0, 24])
    plt.title(' '.join(label)+' '+ ' '.join(choice[0:len(label)])+"\nchannel="+ str(channel) )
    plt.xlabel('time (h)')
    plt.ylabel('fluo (ua)')
    plt.show()
    plt.draw()
    if namePath!='':
        fig1.savefig(namePath)


# In[55]:

label =['CAA1mM']#, 'CAA100uM','CAA400uM', ',CAA700uM', 'CAA0M-eth', 'CAA0M']
high=5
low=0.001
    
#channel='mCherry'
#plotLabelSeveralLog(df, dropMap, ['CAA1mM', 'CAA100uM','CAA0M'],channel, low, high,path+'CAAChannelmCLog.jpg')
#plotLabelSeveralLog(df, dropMap, ['M9-glu1mM', 'M9-glu100uM', 'M9-glu0M'],channel, low, high,path+'M9gluChannelmCLog.jpg')
#plotLabelSeveralLog(df, dropMap, ['pvd'],channel, low, high,path+'pvdChannelmCLog.jpg')

channel='PVD'
#plotLabelSeveralLog(df, dropMap, label,channel, low, high,path+'CAAChannelPVDLog.jpg')
#plotLabelSeveralLog(df, dropMap, ['M9-glu1mM', 'M9-glu100uM', 'M9-glu0M'],channel, low, high,path+'M9gluChannelPVDLog.jpg')
#plotLabelSeveralLog(df, dropMap, ['pvd'],channel, low, high,path+'pvdChannelPVDLog.jpg')

channel='GFP'
plotLabelSeveralLog(df, dropMap, label,channel, low, high,path+'CAAChannelGFPLog.jpg')
#plotLabelSeveralLog(df, dropMap, ['M9-glu1mM', 'M9-glu100uM', 'M9-glu0M'],channel, low, 1,path+'M9gluChannelGFPLog.jpg')
#plotLabelSeveralLog(df, dropMap, ['pvd'],channel, low, high,path+'pvdChannelGFPLog.jpg')


# In[32]:

#plot yield as position

choice=['blue','green','red','cyan','magenta','yellow','black']
colors=choice[0:len(label)];
channel='GFP'
#label=['CAA1mM', 'CAA100uM','CAA0M']
label =['CAA1mM', 'CAA100uM', 'CAA0M','CAA400uM', ',CAA700uM', 'CAA0M-eth', 'CAA0M']


pos=-3
for i,content in enumerate(dropMap[:,1]):
        for j,tmp in enumerate(label):
            if content==tmp:
                p=df[i]
                if channel=='mCherry':
                    plt.scatter(i, y.iloc[pos],c=colors[j], cmap=plt.cm.RdYlGn)
                if channel=='GFP':
                    y=p.fluo_3_mean
                    plt.scatter(i, y.iloc[pos],c=colors[j], cmap=plt.cm.RdYlGn)
                if channel=='PVD':
                    y=p.fluo_1_mean
                    plt.scatter(i, y.iloc[pos],c=colors[j], cmap=plt.cm.RdYlGn)

fig1 = plt.gcf()
#plt.ylim([0, 5])
#plt.xlim([0, 24])
plt.title(' '.join(label)+' '+ ' '.join(choice[0:len(label)])+"\nchannel="+ str(channel) )
plt.xlabel('position')
plt.ylabel('fluo')
plt.show()
plt.draw()
fig1.savefig(path+"yieldPosition"+channel+".png")


# In[16]:


from operator import itemgetter

def movingaverage(interval, window_size):
    window = np.ones(int(window_size))/float(window_size)
    return np.convolve(interval, window, 'same')


def plotLabelSeveralLogAveraged(df, dropMap, label, channel,low=0, high=5, namePath='', avg=20):
    choice=['blue','green','red','cyan','magenta','yellow','black']

    colors=choice[0:len(label)];

    for j,tmp in enumerate(label):
        y=[]
        x=[]
        for i,content in enumerate(dropMap[:,1]):
            if content==tmp:
                    if channel=='mCherry':
                        p=df[i]
                        y.extend(np.log(p.fluo_2_mean))
                        x.extend(p.time/3600)
                    if channel=='GFP':
                        p=df[i]
                        y.extend(np.log(p.fluo_3_mean))
                        x.extend(p.time/3600)
                    if channel=='PVD':
                        p=df[i]
                        y.extend(np.log(p.fluo_1_mean))
                        x.extend(p.time/3600)

        L = sorted(zip(x,y))
        new_x, new_y = zip(*L)
        ysmooth=movingaverage(new_y,avg)
        plt.plot(new_x[avg:],ysmooth[avg:],color=colors[j])

    fig1 = plt.gcf()
        
    #plt.ylim([np.log(low), np.log(high)])
    plt.xlim([0, 30])
    plt.title(' '.join(label)+' '+ ' '.join(choice[0:len(label)])+"\nchannel="+ str(channel) )
    plt.xlabel('time (h)')
    plt.ylabel('fluo (ua)')
    plt.show()
    plt.draw()
    if namePath!='':
        fig1.savefig(namePath)


# In[17]:

label =['CAA1mM', 'CAA100uM','CAA400uM', 'CAA700uM', 'CAA0M-eth', 'CAA0M']
channel='GFP'
plotLabelSeveralLogAveraged(df, dropMap, label,channel, low, high,path+'bactChannelGFPLogAll.jpg',160)
channel='PVD'
plotLabelSeveralLogAveraged(df, dropMap, label,channel, low, high,path+'bactChannelGFPLogAll.jpg',160)


# In[23]:

def func(f,x):
    return f[1]+f[0]*x

def plotDerivative(df, dropMap, label, channel,low=0, high=5, namePath='', nbpt=2):
    choice=['blue','green','red','cyan','magenta','yellow','black']

    colors=choice[0:len(label)];

    for j,tmp in enumerate(label):
        y=[]
        x=[]
        d=[]
        m=0
        R=[]
        
        for i,content in enumerate(dropMap[:,1]):
            if content==tmp:
                    if channel=='mCherry':
                        p=df[i]
                        y=np.array(np.log(p.fluo_2_mean))
                        x=np.array(p.time/3600) 
                    if channel=='GFP':
                        p=df[i]
                        y=np.array(np.log(p.fluo_3_mean))
                        x=np.array(p.time/3600)
                    if channel=='PVD':
                        p=df[i]
                        y=np.array(np.log(p.fluo_1_mean))
                        x=np.array(p.time/3600)
                        
                        
                    #dforth=np.gradient(y[0::nbpt],x[0::nbpt])
                    #dback=np.gradient(y[1::nbpt],x[1::nbpt])   
                    dforth=[]
                    dback=[]
                    for i in range(0,len(x)-1,2): #step two to split back and forth
                        xF=np.array(x[i+0:i+nbpt:2])
                        xB=np.array(x[i+1:i+nbpt+1:2])
                        yF=np.array(y[i+0:i+nbpt:2])
                        yB=np.array(y[i+1:i+nbpt+1:2])

                        
                        f=np.polyfit(xF,yF,1,full=True)
                        dforth.append(f[0][0])
                        #print(f)
                        f=np.polyfit(xB,yB,1,full=True)
                        dback.append(f[0][0])
                        chi_squaredB = np.sum((np.polyval(f[0], xB) - yB) ** 2)/len(xB)
                        R.append(chi_squaredB)
                        chi_squaredF = np.sum((np.polyval(f[0], xF) - yF) ** 2)/len(xF)
                        R.append(chi_squaredF)
                    

                    plt.scatter(x[0:-1:2],np.abs(dforth),color='blue')
                    plt.scatter(x[1::2],np.abs(dback),color='red')
                    plt.plot(x[0:-1:2],np.abs(dforth),color='blue')
                    plt.plot(x[1::2],np.abs(dback),color='red')
                    
                    if m<np.max(dback):
                        m=np.max(dback)
                    if m<np.max(dforth):
                        m=np.max(dforth)
                        

    plt.plot(x,np.ones(len(x))*m*.5)    
    fig1 = plt.gcf()    
    plt.ylim([0, 2])
    plt.xlim([0, 40])
    plt.title(' '.join(label)+' '+ ' '.join(choice[0:len(label)])+"\nchannel="+ str(channel)+"\nnbpts="+str(nbpt) )
    plt.xlabel('time (h)')
    plt.ylabel('derive fluo (ua)')
    plt.show()
    if namePath!='':
        fig1.savefig(namePath)
    
    hist, bins = np.histogram(R, bins=50)
    width = 0.7 * (bins[1] - bins[0])
    center = (bins[:-1] + bins[1:]) / 2
    plt.bar(center, hist, align='center', width=width)
#    plt.xlim([0, 1])
    plt.ylim([0, 10])

    plt.show()



# In[24]:

#plot the derivative of the curves
import warnings
warnings.simplefilter('ignore', np.RankWarning)

label =['CAA1mM']#, 'CAA100uM', 'CAA0M','CAA400uM', ',CAA700uM', 'CAA0M-eth', 'CAA0M']
high=5
low=0.001
for i in range(1,15,2):
    plotDerivative(df, dropMap, label,'GFP', low, high,path+'bactChannelGFPderivLog.jpg',i)


# In[108]:

def func(f,x):
    return f[1]+f[0]*x


def plotDerivativeSeparateBF(df, dropMap, label, channel,low=0, high=5, namePath='', nbpt=2):
    choice=['blue','green','red','cyan','magenta','yellow','black']

    colors=choice[0:len(label)];

    for j,tmp in enumerate(label):
        y=[]
        x=[]
        d=[]
        m=0

        
        for i,content in enumerate(dropMap[:,1]):
            if content==tmp:
                    if channel=='mCherry':
                        p=df[i]
                        y=np.array(np.log(p.fluo_2_mean))
                        x=np.array(p.time/3600) 
                    if channel=='GFP':
                        p=df[i]
                        y=np.array(np.log(p.fluo_3_mean))
                        x=np.array(p.time/3600)
                    if channel=='PVD':
                        p=df[i]
                        y=np.array(np.log(p.fluo_1_mean))
                        x=np.array(p.time/3600)
                        
                        
                    #dforth=np.gradient(y[0::nbpt],x[0::nbpt])
                    #dback=np.gradient(y[1::nbpt],x[1::nbpt])   
                    dforth=[]
                    dback=[]
                    xForth=[]
                    xBack=[]
                    RB=[]
                    RF=[]
                    
                      
                    fig, ax1 = plt.subplots()#s(1,3,sharex=True)
                    fig, ax11 = plt.subplots()#s(1,3,sharex=True)
                    a=np.arange(20)

                    for i in range(0,len(x)-1,2): #step two to split back and forth
                        xF=np.array(x[i+0:i+nbpt:2])
                        xB=np.array(x[i+1:i+nbpt+1:2])
                        yF=np.array(y[i+0:i+nbpt:2])
                        yB=np.array(y[i+1:i+nbpt+1:2])

                        
                        
                        f=np.polyfit(xF,yF,1,full=True)
                        dforth.append(f[0][0])
                        xForth.append(np.mean(xF))
                        
                        thplot=0
                        if xB[0]>thplot and xB[0]<thplot+20: #plot in the exp phase
                            ax1.plot(x,func(f[0],x))
                            plt.ylim(-6, 2)
                            
                        chi_squaredF = np.sum((np.polyval(f[0], xF) - yF) ** 2)/len(xF)
                        RF.append(chi_squaredF)
                        
                        f=np.polyfit(xB,yB,1,full=True)
                        dback.append(f[0][0])
                        xBack.append(np.mean(xB))
                        if xB[0]>thplot and xB[0]<thplot+20: #plot in the exp phase
                            ax11.plot(x,func(f[0],x))
                            plt.ylim(-6, 2)

                        chi_squaredB = np.sum((np.polyval(f[0], xB) - yB) ** 2)/len(xB)
                        RB.append(chi_squaredB)
                        


                    #ax1.scatter(x[0::2],y[0::2],color='blue')#,marker='o', linestyle='-',color='blue')
                    #ax1.scatter(x[1::2],y[1::2],color='red')#,marker='o', linestyle='-',color='red')
                    #ax1.set_ylim(-6, 2)
                    #a=np.arange(20)
                    #ax1.set_title(str(a[:nbpt:2])+' fit')
                    #ax1[0].set(adjustable='box-forced', aspect='equal')
                    #plt.show()

                    #fig, ax1 = plt.subplots()
                    ax1.scatter(x[0::2],y[0::2],color='blue')#,marker='o', linestyle='-',color='blue')
                    ax1.set_ylim([-6, 2])
                    ax2 = ax1.twinx()
                    ax2.plot(xForth,dforth,color='blue',marker='*', linestyle='-')
            
                    ax11.scatter(x[1::2],y[1::2],color='red')#,marker='o', linestyle='-',color='red')
                    ax11.set_ylim([-6, 2])
                    ax2 = ax11.twinx()
                    ax2.plot(xBack,dback,color='red',marker='*', linestyle='-')
                    
                    
                    #ax1[1].set(adjustable='box-forced', aspect='equal')
                    #ax2.set(adjustable='box-forced', aspect='equal')
                    ax1.set_title(str(a[:nbpt:2])+' derivative')
                    ax11.set_title(str(a[:nbpt:2])+' derivative')
                    #plt.show()
                    

                    #ax1.plot(x[0::2],y[0::2],marker='o', linestyle='-',color='blue')
                    #ax3 = ax1.twinx()
                    #ax3.plot(x[0:-1:2],RB,marker='.', linestyle='--',color='blue')
                    
                    #ax11.plot(x[1::2],y[1::2],marker='o', linestyle='-',color='red')
                    #ax3 = ax11.twinx()
                    #ax3.plot(x[1::2],RF,marker='.', linestyle='--',color='red')
                    
                    #ax1[2].set(adjustable='box-forced', aspect='equal')
                    #ax3.set(adjustable='box-forced', aspect='equal')
                    
                    #a=np.arange(20)
                    #ax1.set_title('residual '+str(a[:nbpt:2])+str())
                    #fig.tight_layout()
                    #ax1.set_ylim([-6, 2])
                    #ax3.set_ylim([0, .5])
     
                    plt.show()

                    if m<np.max(dback):
                        m=np.max(dback)
                    if m<np.max(dforth):
                        m=np.max(dforth)
                        

    #plt.plot(x,np.ones(len(x))*m*.5)  
    #fig1 = plt.gcf()    
    #plt.ylim([0, 2])
    #plt.xlim([0, 40])
    #plt.title(' '.join(label)+' '+ ' '.join(choice[0:len(label)])+"\nchannel="+ str(channel)+"\nnbpts="+str(nbpt) )
    #plt.xlabel('time (h)')
    #plt.ylabel('derive fluo (ua)')
    #plt.show()
    #if namePath!='':
        #fig1.savefig(namePath)
    
    #hist, bins = np.histogram(R, bins=50)
    #width = 0.7 * (bins[1] - bins[0])
    #center = (bins[:-1] + bins[1:]) / 2
    #plt.bar(center, hist, align='center', width=width)
#    plt.xlim([0, 1])
    #plt.ylim([0, 10])

    #plt.show()



# In[109]:

import warnings
warnings.simplefilter('ignore', np.RankWarning)
print('start')
label =['CAA1mM']#, 'CAA100uM', 'CAA0M','CAA400uM', ',CAA700uM', 'CAA0M-eth', 'CAA0M']
high=5
low=0.001
plotDerivativeSeparateBF(df, dropMap, label,'GFP', low, high,path+'bactChannelGFPderivLog.jpg',7)


# In[145]:

def plotDerivativeSeparate(df, dropMap, label, channel,low=0, high=5, namePath='', nbpt=2):
    choice=['blue','green','red','cyan','magenta','yellow','black']

    colors=choice[0:len(label)];

    for j,tmp in enumerate(label):
        y=[]
        x=[]
        d=[]
        m=0

        
        for i,content in enumerate(dropMap[:,1]):
            if content==tmp:
                    if channel=='mCherry':
                        p=df[i]
                        y=np.array(np.log(p.fluo_2_mean))
                        x=np.array(p.time/3600) 
                    if channel=='GFP':
                        p=df[i]
                        y=np.array(np.log(p.fluo_3_mean))
                        x=np.array(p.time/3600)
                    if channel=='PVD':
                        p=df[i]
                        y=np.array(np.log(p.fluo_1_mean))
                        x=np.array(p.time/3600)
                        
                        
                    #dforth=np.gradient(y[0::nbpt],x[0::nbpt])
                    #dback=np.gradient(y[1::nbpt],x[1::nbpt])   
                    dforth=[]
                    dback=[]
                    xForth=[]
                    xBack=[]
                    xBF=[]
                    dBF=[]
                    RB=[]
                    RF=[]
                    
                      
                    fig, ax1 = plt.subplots()#s(1,3,sharex=True)
                    a=np.arange(20)

                    for i in range(0,len(x)-1,2): #step two to split back and forth
                        xF=np.array(x[i+0:i+nbpt:2])
                        xB=np.array(x[i+1:i+nbpt+1:2])
                        yF=np.array(y[i+0:i+nbpt:2])
                        yB=np.array(y[i+1:i+nbpt+1:2])

                        
                        
                        f=np.polyfit(xF,yF,1,full=True)
                        dforth.append(f[0][0])
                        xForth.append(np.mean(xF))
                        
                        thplot=0
                        if xB[0]>thplot and xB[0]<thplot+20: #plot in the exp phase
                            ax1.plot(x,func(f[0],x))
                            plt.ylim(-6, 2)
                            
                        chi_squaredF = np.sum((np.polyval(f[0], xF) - yF) ** 2)/len(xF)
                        RF.append(chi_squaredF)
                        
                        f=np.polyfit(xB,yB,1,full=True)
                        dback.append(f[0][0])
                        xBack.append(np.mean(xB))
                        if xB[0]>thplot and xB[0]<thplot+20: #plot in the exp phase
                            ax1.plot(x,func(f[0],x))
                            plt.ylim(-6, 2)

                        chi_squaredB = np.sum((np.polyval(f[0], xB) - yB) ** 2)/len(xB)
                        RB.append(chi_squaredB)
                        


                    #ax1.scatter(x[0::2],y[0::2],color='blue')#,marker='o', linestyle='-',color='blue')
                    #ax1.scatter(x[1::2],y[1::2],color='red')#,marker='o', linestyle='-',color='red')
                    #ax1.set_ylim(-6, 2)
                    #a=np.arange(20)
                    #ax1.set_title(str(a[:nbpt:2])+' fit')
                    #ax1[0].set(adjustable='box-forced', aspect='equal')
                    #plt.show()

                    #fig, ax1 = plt.subplots()
                    ax1.scatter(x[0::2],y[0::2],color='blue')#,marker='o', linestyle='-',color='blue')
                    ax1.scatter(x[1::2],y[1::2],color='red')#,marker='o', linestyle='-',color='red')
                    ax1.set_ylim([-6, 2])
                    ax2 = ax1.twinx()
                    xBF=xForth+xBack
                    xBF[::2]=xBack
                    xBF[1::2]=xForth
                    dBF=dforth+dback
                    dBF[::2]=dback
                    dBF[1::2]=dforth
                    ax2.plot(xBF,dBF,color='green',marker='*', linestyle='-')
                    ax2.set_ylim([0, 1])

     
                    plt.show()

                    if m<np.max(dback):
                        m=np.max(dback)
                    if m<np.max(dforth):
                        m=np.max(dforth)
                        

    #plt.plot(x,np.ones(len(x))*m*.5)  
    #fig1 = plt.gcf()    
    #plt.ylim([0, 2])
    #plt.xlim([0, 40])
    #plt.title(' '.join(label)+' '+ ' '.join(choice[0:len(label)])+"\nchannel="+ str(channel)+"\nnbpts="+str(nbpt) )
    #plt.xlabel('time (h)')
    #plt.ylabel('derive fluo (ua)')
    #plt.show()
    #if namePath!='':
        #fig1.savefig(namePath)
    
    #hist, bins = np.histogram(R, bins=50)
    #width = 0.7 * (bins[1] - bins[0])
    #center = (bins[:-1] + bins[1:]) / 2
    #plt.bar(center, hist, align='center', width=width)
#    plt.xlim([0, 1])
    #plt.ylim([0, 10])

    #plt.show()


# In[146]:

import warnings
warnings.simplefilter('ignore', np.RankWarning)
print('start')
label =['CAA1mM']#, 'CAA100uM', 'CAA0M','CAA400uM', ',CAA700uM', 'CAA0M-eth', 'CAA0M']
high=5
low=0.001
plotDerivativeSeparate(df, dropMap, label,'GFP', low, high,path+'bactChannelGFPderivLog.jpg',7)


# In[34]:


from scipy.optimize import curve_fit

def sigmoid(x, k, x0, A, B):
    y =  B+A / (1 + np.exp(-k*(x-x0)))
    return y


def plotSig(df, dropMap, label, channel,low=0, high=5):
    choice=['blue','green','red','cyan','magenta','yellow','black']

    colors=choice[0:len(label)];

    for j,tmp in enumerate(label):
        y=[]
        x=[]
        d=[]
        m=0

        
        for i,content in enumerate(dropMap[:,1]):
            if content==tmp:
                    if channel=='mCherry':
                        p=df[i]
                        y=np.array(np.log(p.fluo_2_mean))
                        x=np.array(p.time/3600) 
                    if channel=='GFP':
                        p=df[i]
                        y=np.array(np.log(p.fluo_3_mean))
                        x=np.array(p.time/3600)
                    if channel=='PVD':
                        p=df[i]
                        y=np.array(np.log(p.fluo_1_mean))
                        x=np.array(p.time/3600)
                    
                    y=y-np.mean(y[0:1])
                    
                    xB=x[0::2]
                    yB=y[0::2]
                    
                    xF=x[1::2]
                    yF=y[1::2]
                    
                    #sigmoid(x, k, x0, A, B)
                    boundDown=[.1,5.,5., -.1]
                    boundUp=[2.,20.,7., .1]

                    poptB, pcovB = curve_fit(sigmoid,xB,yB,bounds=(boundDown, boundUp) )
                    poptF, pcovF = curve_fit(sigmoid, xF, yF,bounds=(boundDown, boundUp))
                    print(poptB, poptF)
                      
                    fig, ax1 = plt.subplots()#s(1,3,sharex=True)
 
                    ax1.scatter(x[0::2],y[0::2],color='blue')#,marker='o', linestyle='-',color='blue')
                    ax1.scatter(x[1::2],y[1::2],color='red')#,marker='o', linestyle='-',color='red')
                    #ax1.set_ylim([-6, 2])
                    
                    ax1.plot(xF, sigmoid(xF, *poptF),color='green', linestyle='-')
                    ax1.plot(xB, sigmoid(xB, *poptB),color='green', linestyle='-')
    

     
                    plt.show()

          


# In[36]:

import warnings
warnings.simplefilter('ignore', np.RankWarning)
print('start')
#label =['CAA1mM']#, 'CAA100uM', 'CAA0M','CAA400uM', ',CAA700uM', 'CAA0M-eth', 'CAA0M']
label = ['CAA-1']#, 'CAA-5', 'CAA-10', 'CAA-25', 'CAA-50', 'CAA']
high=5
low=0.001
plotSig(df, dropMap, label,'GFP', low, high)


# In[7]:

get_ipython().system('jupyter nbconvert --to script dropAnalysis.ipynb ')
get_ipython().system('jupyter rename dropAnalysis.ipynb path+dropAnalysis.ipynb')

