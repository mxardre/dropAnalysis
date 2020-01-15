#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 20 21:31:36 2019

@author: admin
"""

import pandas as pd
import numpy as np
import os
import shutil
from os import listdir
from os.path import isfile, join, isdir
import matplotlib.pyplot as plt
import matplotlib
import scipy.optimize
import scipy.interpolate as si
#import seaborn as sns

#from ggplot import *

from pandas import DataFrame as pddf
from pandas import Series as pds

import bisect
import csv
import glob

import sys

module_path = os.path.abspath(os.path.join('/Users/bacterie/Documents/Programme/fitderivpackage1.02'))
if module_path not in sys.path:
    sys.path.append(module_path)
from fitderiv import fitderiv 


def getValueInLabel(label,path):

    
    p=path.split('/')
    posOfFolderInSource=-4
#try:
    if 'Exp_Incub_4h30-6h30-24h-52h' == p[posOfFolderInSource] :
        #c = label.split('SBW25GFP-')
        #d = c[1]
        #b = d.split('h N=')
        c=label.split('h')
        b=c[0]

    if 'Exp_Dilution_PvdS' == p[posOfFolderInSource] :
        c = label.split('CAA-PvdS-')
        d = c[1]
        e = d.split(' N')
        b=e[0]

    if 'Exp_DilutionWTstock2' == p[posOfFolderInSource] :
        c = label.split('CAA-WT-')
        e = c[1]
        d = e.split(' N')  
        b=d[0]

    if 'Exp_Dilution_5_50_500_3500' == p[posOfFolderInSource] :
        c = label.split('CAA-WT-')
        d = c[1]
        e = d.split(' N')
        b=e[0]

    if 'Exp_IronCAA_WT' == p[posOfFolderInSource] :
        c=label.split('CAA-WT-')
        d=c[1]
        e = d.split('uM')
        b=e[0]

    if 'Exp_BipyCAA_WT' == p[posOfFolderInSource] :
        c=label.split('uM')
        b=c[0]

    if 'Exp_BipyCAA_PvdS' == p[posOfFolderInSource] :
        c=label.split('uM')
        b=c[0]

    if 'Exp_Incubation_14h-19h-43h' == p[posOfFolderInSource] :
        c=label.split('SBW25GFP-')
        d=c[1]
        e = d.split('h')
        b=e[0]

    if 'Exp_Incub_4h30-6h30-24h-52h' == p[posOfFolderInSource] :
        c=label.split('SBW25-WT_CAA_')
        d=c[1]
        e = d.split('h')
        b=e[0]
        
    if 'Exp_Incubation_4h30_6h30_24h_43h13_1bactPerDrop' == p[posOfFolderInSource] :
        c=label.split('SBW25-WT_CAA_')
        d=c[1]
        if d[:3] == 'pur':
            b=0
        else :
            e = d.split('h')
            b=e[0]

        
    if 'Exp_Incubation_merge' == p[posOfFolderInSource] :
        f=p[posOfFolderInSource+1] 
        ff=f.split('dilution_CAA_WT-')
        
        if ff[-1]== 'stock24h30-6h30-24h-52h':
            c=label.split('SBW25-WT_CAA_')
            d=c[1]
            e = d.split('h')
            b=e[0]
            
        if ff[-1]== 'stock14h33-19h-43h':
            c=label.split('SBW25GFP-')
            d=c[1]
            e = d.split('h')
            b=e[0]

    if 'Exp_IronPvd' == p[posOfFolderInSource] :
        
        c=label.split('CAA-pur-')
        d=c[1]
        e = d.split('uM')
        b=e[0]
        
    if 'Exp_IronCAA_PvdS' == p[posOfFolderInSource] :
        c=label.split('CAA-PvdS-')
        d=c[1]
        e = d.split('uM')
        b=e[0]
        
    if 'Exp_Dilution_M9glycerol' == p[posOfFolderInSource] :
        c=label.split('SBW25-WT_M9gly_')
        b=c[1]
        
    if 'Exp_GreenRed' == p[posOfFolderInSource] :
        val=label

    else:
        val= float(b)
        
    return val
#except:

    #print('getValueInLabel() label not found')
    #print(path)
    #print(label)
    #print( p[posOfFolderInSource])




def findIdx(your_list,item):
    lst=pddf(your_list)
    lst2=lst[0].values
    lst3=lst2.tolist()
    return lst3.index(item)

def getOrderLabel(orderList,toOrganize):
    for idx in toOrganize:
        order.append(orderList.index(idx))
    return order

def organizeLabel(order,toOrganize):
    zipped=zip(toOrganize,order)
    organized=sorted(zipped, key=lambda x: x[1])
    l=zip(*organized)
    return list(l)[0]




def loadData(path):
    drpfiles=sorted([path+'droplets/' + f for f in listdir(path + 'droplets/') if isfile(join(path + 'droplets/',f)) and f!='.DS_Store'])
    i=0
    df={}
    for file in drpfiles:
        df[i]=pd.read_csv(file)
        i+=1 

    tpfile=pd.read_csv(path+'droplet.csv')
    dropMap=[]
    label=list(set(tpfile['group']))
    label2=[l for l in label if l!='Empty' and l!='CAA']
    dropMap=np.array(list(zip(tpfile['well'],tpfile['group'])))
    
    return [dropMap, df, label2]




def getfitData(df, j, channel,startTime,timeThresh,threshOutlayer, incCarte): 

    y=[]
    x=[] 
    seuilDetectionGFP=0#2e-2V checked on an empty droplet
    #seuilDetectionPVD=3e-2V checked on an empty droplet codé en dure
    
    
    try:
        if channel=='SpeedSize':
    
            p=df[j]
            idxThresDetect=1 #find threshold of data above detection 1e-2V
            y=np.array(p.speed/p.size)
            x=np.array(p.time/3600) 
            x=x-x[0] #reset time
            ystd=[]
            
            #remove nan value            
            idxNan=[i for i, j in enumerate(np.isnan(y)) if j] 
            x = np.delete(x,idxNan)
            y = np.delete(y,idxNan)
            ystd = np.delete(ystd,idxNan) 
        
        if channel=='RFP':
    
            p=df[j]
            idxThresDetect=bisect.bisect(p.fluo_2_median,seuilDetectionGFP)#find threshold of data above detection 5e-3V
            y=np.array( np.log(p['fluo_2_area']*p['speed']/p['size']) - np.log(p['fluo_2_area'][1]*p['speed'][1]/p['size'][1]) )
            x=np.array(p.time/3600) 
            x=x-x[0] #reset time
            ystd=np.array(np.divide(p['fluo_2_std']+incCarte/2,p['fluo_2_median']))
    
            #remove nan value            
            idxNan=[i for i, j in enumerate(np.isnan(y)) if j] 
            x = np.delete(x,idxNan)
            y = np.delete(y,idxNan)
            ystd = np.delete(ystd,idxNan) 
            
        if channel=='GFP':
    
            
            p=df[j]
            idxThresDetect=bisect.bisect(p.fluo_3_median,seuilDetectionGFP)#find threshold of data above detection 5e-3V
            y=np.array( np.log(p['fluo_3_area']*p['speed']/p['size']) - np.log(p['fluo_3_area'][1]*p['speed'][1]/p['size'][1]) )
            
            
            x=np.array(p.time/3600,dtype=np.float64)
            x=x-x[0] #reset time
            ystd=np.array(np.divide(p['fluo_3_std']+incCarte/2,p['fluo_3_median']))
            
            #remove nan value            
            idxNan=[i for i, j in enumerate(np.isnan(y)) if j] 
            x = np.delete(x,idxNan)
            y = np.delete(y,idxNan)
            ystd = np.delete(ystd,idxNan)    
            
        if channel=='PVD':
    
            p=df[j]
            idxThresDetect=bisect.bisect(p.fluo_1_median,seuilDetectionGFP)#find threshold of data above detection 6e-3V
            y=np.array( np.log(p['fluo_3_area']*p['speed']/p['size']) - np.log(p['fluo_1_area'][1]*p['speed'][1]/p['size'][1]) )
            x=np.array(p.time/3600,dtype=np.float64)
            x=x-x[0] #reset time
            ystd=np.array(np.divide(p['fluo_1_std']+incCarte/2,p['fluo_1_median']))

            #remove nan value            
            idxNan=[i for i, j in enumerate(np.isnan(y)) if j] 
            x = np.delete(x,idxNan)
            y = np.delete(y,idxNan)
            ystd = np.delete(ystd,idxNan) 
            
        if channel=='PVDoverGFP':
            
            p=df[j]
            
            idxThresDetectGFP=bisect.bisect(p.fluo_3_median,seuilDetectionGFP)#find threshold of data above detection 5e-3V
            idxThresDetectPVD=bisect.bisect(p.fluo_1_median,seuilDetectionGFP)
            
            idxThresDetect= max(idxThresDetectPVD,idxThresDetectGFP)
            
            yFluo=np.array(np.log(p['fluo_3_area']*p['speed']/p['size']),dtype=np.float64) #to screen for outlayer
            y=np.array(p['fluo_1_area']/p['fluo_3_area']) #no need of log here
            x=np.array(p.time/3600,dtype=np.float64)
            x=x-x[0] #reset time
            ystd=np.array(np.divide(p['fluo_1_std']+incCarte/2,p['fluo_1_median']),dtype=np.float64)
            
            #remove nan value            
            idxNan=[i for i, j in enumerate(np.isnan(y)) if j] 
            x = np.delete(x,idxNan)
            y = np.delete(y,idxNan)
            ystd = np.delete(ystd,idxNan) 

        #filtre outlayer
        #the outlayer are detected with the median of the signal
        #the threshold of detection is detected with the median of the signal. The minimum detection is 2e-2V 
        #signal above empty droplets. The plotted signal is an integrated value of the fluo (it is not comparable to 2e-2V directly)
        
        idxStart=max(idxThresDetect,bisect.bisect(x, startTime))
    
        idxThres=bisect.bisect(x, min(max(x),timeThresh))
        idxOutlayer=bisect.bisect(x, min(max(x),threshOutlayer))
        
        if len(y)>idxOutlayer:
            if channel=='GFP'or channel=='PVD' or channel=='RFP':
                outlayer= y[idxOutlayer]>y[2]+1 or threshOutlayer<=0
               
            elif channel=='PVDoverGFP':
                outlayer= yFluo[idxOutlayer]>yFluo[2]+2 or threshOutlayer<=0 
            elif channel=='SpeedSize':
                outlayer= True
            else :
                print('NO CHANNEL FOUND')
    
        else:
            print('BUG idxOutlayer not long enough')
            outlayer=False
    
                
    except Exception as inst: 
         outlayer = False
         exc_type, exc_obj, exc_tb = sys.exc_info()
         fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
         print(exc_type, fname, exc_tb.tb_lineno)
         print('bug checkdata')
         print(type(inst))     # the exception instance
         print(inst.args)     # arguments stored in .args        


    if outlayer :
        return [x[idxStart:idxThres], y[idxStart:idxThres],True]
    else:
        print('outlayer')
        return [[],[],False]




def checkfitData(df, dropMap, label, channel,folder,startTime,timeThresh,threshOutlayer, incCarte): 
    for k,dataFile in enumerate(label):
        

        fig=plt.figure()
        date=folder.split('/')
        plt.title(date[-3]+' '+dataFile)
        #plt.ylim(np.log([1e-4,5]))
        plt.xlim([0,timeThresh])
        
        err=np.array([])
        for j,content in enumerate(dropMap[:,1]):
            y=[]
            x=[] 


            if content==dataFile: 
                
                [x,y,flag] = getfitData(df, j, channel,startTime,timeThresh,threshOutlayer, incCarte) 
                
                if flag==True:
                    if channel=='RFP':
                        c='r'            
                    if channel=='GFP':
                        c='b'                    
                    if channel=='PVD':
                        c='g'
                    if channel=='PVDoverGFP':
                        c='k'
                        plt.ylim([0,3])
                    if channel=='SpeedSize':   
                        c='m'
                        a=100*np.std(y)/np.mean(y)
                        err=np.concatenate([err,[a]])
                        plt.title(date[-3]+' '+dataFile+' std/m='+str("{0:.2g}".format(np.mean(err)))+'%')
                        plt.ylim([np.min(y)/2,2*np.max(y)])
                        
                    
                    plt.scatter(x, y,color=c)
                #else: 
                    #print("data removed by thresh drop="+str(j))

        plt.grid()
        plt.show()  
        print(np.mean(np.mean(err)))
        fig.savefig(folder+dataFile+'checkDataFit'+channel+'.jpg')




def fitDataIndiv(df, dropMap, label, channel,folder,startTime,timeThresh,threshOutlayer, incCarte, display=0, deletionData=False): 

    print('RESET TIME!')
    pathResults = folder+'resultIndiv/'
    if not os.path.exists(pathResults):
        os.makedirs(pathResults)
    else :
        if deletionData==True:
            shutil.rmtree(pathResults)
            os.makedirs(pathResults)


    
    for j,dataFile in enumerate(label):
        
        pathFitPlot=folder+'resultIndiv/'+dataFile+'/'
        if not os.path.exists(pathFitPlot):
            os.makedirs(pathFitPlot)
        else :
            if deletionData==True:
                shutil.rmtree(pathResults)
                os.makedirs(pathResults)
                
        for j,content in enumerate(dropMap[:,1]):
            y=[]
            x=[]

            if content==dataFile:    
                if channel=='RFP':
                    c='r'            
                if channel=='GFP':
                    c='b'                    
                if channel=='PVD':
                    c='g'
                if channel=='PVDoverGFP':
                    c='b'
                if channel=='SpeedSize':   
                    c='m'
                
                [x,y,flag] = getfitData(df, j, channel,startTime,timeThresh,threshOutlayer, incCarte) 
                
                
                
                if flag==True:
                    idx=np.where(np.isinf(y))
                    y=np.delete(y,idx[0])
                    x=np.delete(x,idx[0])
                    
                    try:
                        # redefine bounds and run inference
                        #hyperparameters 
                        #0: amplitude
                        #1: flexibility
                        #2: measurement error
                        #initial    b= {0: [-1,5], 1: [-5,2], 2: [-7,1]}
                        b= {0: [-1,5], 1: [-5,0], 2: [-10,0]}

                        #apply on the range from start to max time
                        q= fitderiv(x, y, bd= b, logs= False,esterrs=False)

                        print(q.ds['max df'])
                        if display==True:
                            
                            # plot results
                            fig=plt.figure()
                            #plt.subplot(2,1,1)
                            q.plotfit('f')
                            plt.grid()
                            plt.title(dataFile)
                            plt.show()
                            fig.savefig(pathFitPlot+dataFile+'fitIndiv_Drop'+str(j)+channel+'.jpg')
                            #plt.subplot(2,1,2)
                            fig=plt.figure()
                            q.plotfit('df')
                            plt.grid()
                            plt.title(dataFile)
                            plt.show()
                            fig.savefig(pathFitPlot+dataFile+'derivativeIndiv_Drop'+str(j)+channel+'.jpg')


                        # export results
                        statdict=q.printstats()
                        with open(folder+'resultIndiv/'+dataFile+'resultIndiv_Drop'+str(j)+channel+'.csv','w') as csv_file:
                            writer = csv.writer(csv_file)
                            for key, value in statdict.items():
                                writer.writerow([key, value])


                    except Exception as inst:
                        print('bug od')
                        print(type(inst))     # the exception instance
                        print(inst.args)     # arguments stored in .args




def barPlot(labelList,y,stdy,ymin,ymax,title,ylabel,path):
    # Build the plot
    fig, ax = plt.subplots()
    x_pos = np.arange(len(labelList))
    ax.bar(x_pos, y, yerr=stdy, align='center', alpha=0.5, ecolor='black', capsize=10)
    ax.set_ylabel(ylabel)
    ax.set_xticks(x_pos)
    plt.xticks(rotation=90)
    ax.set_xticklabels(labelList)
    ax.set_ylim([ymin,ymax])
    ax.set_title(title)
    ax.yaxis.grid(True)

    # Save the figure and show
    plt.tight_layout()
    plt.savefig(path,bbox_inches='tight', format='pdf')
    plt.show()
    
def boxPlot(labelList,y,ymin,ymax,title,ylabel,path):
    # Build the plot
    fig, ax = plt.subplots()
    x_pos = np.arange(len(labelList))
    plt.boxplot(y)
    ax.set_ylabel(ylabel)
    ax.set_xticks(x_pos)
    plt.xticks(rotation=90)
    ax.set_xticklabels(labelList)
    ax.set_ylim([ymin,ymax])
    ax.set_title(title)
    ax.yaxis.grid(True)

    # Save the figure and show
    plt.tight_layout()
    plt.savefig(path,bbox_inches='tight', format='pdf')
    plt.show()




def findIdx(your_list,item):

    lst=pddf(your_list)
    lst2=lst[0].values
    lst3=lst2.tolist()
    return lst3.index(item)




def getDataHistogram(label,path,channel):

    folder=path+'resultIndiv/'
    [dropMap, df, nn]=loadData(path)

    stdgRate = pddf()
    gRate = pddf()
    lag = pddf()
    stdlag = pddf()
    yld = pddf()
    stdyield = pddf()
    
    for n,labelName in enumerate(label):

        stdgRateList = []
        gRateList = []
        lagList = []
        stdlagList = []
        yieldList = []
        stdyieldList = []
 
        for file in os.listdir(folder):
            if file.startswith(labelName+'resultIndiv_'):
                if file.endswith(channel+".csv"):

                    try:  
                        with open(folder+file, 'r') as f:
                            reader = csv.reader(f)
                            your_list = list(reader)

                        stdgRateList.append(float(your_list[findIdx(your_list,'max df std')][1]))
                        gRateList.append(float(your_list[findIdx(your_list,'max df')][1]))
                        lagList.append(float(your_list[findIdx(your_list,'lag time')][1]))
                        stdlagList.append(float(your_list[findIdx(your_list,'lag time std')][1]))
                        yieldList.append(float(your_list[findIdx(your_list,'max y')][1]))
                        stdyieldList.append(float(your_list[findIdx(your_list,'max y std')][1]))

                    except Exception as inst:
                        donothing=0
                        print(labelName)
                        print(folder+file)
                        print(type(inst))     # the exception instance
        
        df = pddf({labelName:gRateList})
        gRate = pd.concat([gRate,df], axis=1)
        
        df = pddf({labelName:lagList})
        lag = pd.concat([lag,df], axis=1)
             
        df = pddf({labelName:yieldList})
        yld = pd.concat([yld,df], axis=1) 
        
        df = pddf({labelName:stdgRateList})
        stdgRate = pd.concat([stdgRate,df], axis=1)
              
        df = pddf({labelName:stdlagList})
        stdlag = pd.concat([stdlag,df], axis=1)
             
        df = pddf({labelName:stdyieldList})
        stdyield = pd.concat([stdyield,df], axis=1)
        
    
    return [stdgRate, gRate, lag, stdlag, yld, stdyield]


def poolDataInterpolate(dropMap, df, label,channel, path, incCarte=5e-3, nbReps=80, startTime=2, timeThresh=40, threshOutlayer=17, display= False): 
    for j,tmp in enumerate(label):
        dataX=np.array([])
        dataY=np.array([])
        data=pddf()

        print(tmp)
        nbAdded=0
        p=df[0]
        x0=np.array(p.time/3600,dtype=np.float64)

        p=df[len(df)-4]#antepenultien run pour creer le range de temps sur lequel interpolle
        xLast=np.array(p.time/3600,dtype=np.float64)
        #x2=np.arange(.5, min(max(x0)-min(x0),max(xLast)-min(xLast)), .5)
        x2=np.arange(.5, timeThresh, .5)
        data['time']=pd.Series(x2).values

        for j,content in enumerate(dropMap[:,1]):

            y=[]
            x=[]

            if content==tmp:    
                if channel=='RFP':
                    p=df[j]
                    y=np.array(np.log(p.fluo_2_area*p.speed/p.size)-np.log(p.fluo_2_area[0]*p.speed[0]/p.size[0]))
                    x=np.array(p.time/3600) 
                    x=x-x[0] #reset time
                    ystd=np.array(np.divide(p.fluo_2_std+incCarte/2,p.fluo_2_median))
                if channel=='GFP':
                    p=df[j]
                    y=np.array(np.log(p.fluo_3_area*p.speed/p.size)-np.log(p.fluo_3_area[0]*p.speed[0]/p.size[0]))
                    x=np.array(p.time/3600,dtype=np.float64)
                    x=x-x[0] #reset time
                    ystd=np.array(np.divide(p.fluo_3_std+incCarte/2,p.fluo_3_median))

                if channel=='PVD':
                    p=df[j]
                    y=np.array(np.log(p.fluo_1_area*p.speed/p.size)-np.log(p.fluo_1_area[0]*p.speed[0]/p.size[0]))
                    x=np.array(p.time/3600,dtype=np.float64)
                    x=x-x[0] #reset time
                    ystd=np.array(np.divide(p.fluo_1_std+incCarte/2,p.fluo_1_median))

                if channel=='PVDoverGFP':
                    p=df[j]
                    y=np.array(p.fluo_1_area/(p.fluo_3_area)) #no need of log here
                    x=np.array(p.time/3600,dtype=np.float64)
                    x=x-x[0] #reset time
                    ystd=np.array(np.divide(p.fluo_1_std+incCarte/2,p.fluo_1_median),dtype=np.float64)


                #filtre outlayer
                #idxStart=bisect.bisect(x, startTime)
                #idxThres=bisect.bisect(x, min(max(x),timeThresh))
                idxOutlayer=bisect.bisect(x, threshOutlayer)

                if len(y)>idxOutlayer:
                    if y[idxOutlayer]>y[2]+2 or threshOutlayer<=0 :
                        # measure at run 0 does not count
                        f = si.interp1d(x, y)
                        try:

                            y2=f(x2)
                            data[str(j)]=pd.Series(y2).values
                            nbAdded+=1

                        except Exception as inst:
                            print('poolData() remove this data')
                            print(j)
                            print(type(inst))     # the exception instance
                            print(inst.args)     # arguments stored in .args
                            print(inst)
                            print('x2')
                            print(x2)
                            print('y')
                            print(y)
                            print('x')
                            print(x)

                #plot the growth curve with the fit and the derivative
                if display == True:
                    pathgrowth = path+'growth/'
                    if not os.path.exists(pathgrowth):
                        os.makedirs(pathgrowth)

                    #plt.plot(x2,y2,marker='o', linestyle='-',color='blue')
                    plt.plot(x,y,marker='o', linestyle='-',color='blue')

                    #ax1.set_ylim([np.log(low), np.log(high)])
                    #ax1.set_yticks(range(np.int(np.log(low)), np.int(np.log(high))))
                    plt.ylabel('log(fluo)')
                    plt.xlabel('time (h)')

                    plt.title(tmp +' drp'+str(j)+' area')
                    plt.savefig(pathgrowth+'drp'+str(j)+'_deriveAll_nbpt'+'_'+tmp, format='pdf')
                    plt.show()

                if nbAdded>nbReps:
                    break
            
        data.to_csv(path+tmp+channel+'Interp.csv', index=False )


        


def checkDataInterpolated(folder,label,channel):
    for dataFile in label:
        file=folder+dataFile+channel+'Interp.csv'
        data=pd.read_csv(file)
        # load data
        t=np.array(data['time'])
        fig=plt.figure()

        if channel=='GFP':
            c='b'
        if channel=='RFP':
            c='r'
        if channel=='PVD':
            c='g'
        if channel=='PVDoverGFP':
            c='k'

        
        for column in data:
            if column!='time':
                od= np.array(data[column])
                plt.scatter(t,od, color=c)
                plt.title(dataFile)
#        plt.ylim([-10,-2])
        plt.show()
        fig.savefig(folder+dataFile+'checkDataInterp'+channel+'.jpg',bbox_inches='tight')
    
def plot_distribution(ax,key,base,val,m_path):
    
    if key==5 or key==4.5 or key==4 or key==14:
        if '2018-02-21_dilution_CAA_WT-stock2' in m_path:
            ax.plot(base[:-1], val, '-*k', markerSize=10)
        else:
            ax.plot(base[:-1], val, '-k')
    elif key==50 or key==25 or key==6.5 or key==6:
        ax.plot(base[:-1], val, '-b')
    elif key==500 or key==100 or key==24  or key==19:
        ax.plot(base[:-1], val, '-r')
    elif key==3500 or key==1000 or key==52 or key==43:
        ax.plot(base[:-1], val, '-g')                                    
    else :
        print('not printed:'+str(key))
        
def my_legend(rootPath):
        
        if 'Exp_Incubation_4h30_6h30_24h_43h13_1bactPerDrop' in rootPath:
            lgd= 'black=4h30 \nblue=6h30 \nred=24h \ngreen=43h'
            
        elif 'Exp_Incub_4h30-6h30-24h-52h' in rootPath:
            lgd='black=4h30 \nblue=6h30 \nred=24h \ngreen=52h'
            
        elif 'Exp_DilutionWTstock2' in rootPath:
            lgd='blue=25 \nred=100 \ngreen=1000'
            
        elif 'Exp_Dilution_5_50_500_3500' in rootPath:
            lgd='black=5 \nblue=50 \nred=500 \ngreen=3500'
            
        elif 'Exp_Incubation_14h-19h-43h' in rootPath:
            lgd='black=14h \nred=19h \ngreen=43h'
            
        elif 'Exp_Dilution_M9glycerol' in rootPath:
            lgd='black=1 \nblue=10 \nred=100 \ngreen=1000'
            
        elif 'Exp_Dilution_PvdS' in rootPath:
            lgd='black=5 \nblue=50 \nred=500 \ngreen=3500'
        
        else :
            lgd='Nan'
            
        return lgd
