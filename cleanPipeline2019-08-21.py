#!/usr/bin/env python
# coding: utf-8

# In[55]:

import os
import sys

module_path = os.path.abspath(os.path.join('/Users/bacterie/Documents/Experiences/Milli/jupyterScript/'))
if module_path not in sys.path:
    sys.path.append(module_path)
from functionsCleanPipeline import *




#define the folder in which the good experiement are
source='/Users/bacterie/Documents/Experiences/Milli/'
#rootPath='/Users/maxime/Documents/experiences/milli/expDilution_5_50_500_3500/2018-02-09_dilution_CAA_WT-stock2/'
#folder=rootPath+'analysis/'
#rootPath='/Users/maxime/Documents/experiences/milli/goodExperiments/'
#rootPath='/Users/maxime/Documents/experiences/milli/experienceDilutionWTstock2/'


#rootPath='/Users/maxime/Documents/experiences/milli/expDilution_5_50_500_3500/'

#rootPath='/Users/maxime/Documents/experiences/milli/experienceGreenRed/'
#rootPath='/Users/maxime/Documents/experiences/milli/expBipyCAA_WT/'
#rootPath='/Users/maxime/Documents/experiences/milli/expBipyCAA_PvdS/'

#rootPathList=['/Users/maxime/Documents/experiences/milli/expBipyCAA_WT/', '/Users/maxime/Documents/experiences/milli/expBipyCAA_PvdS/']

#rootPath=source+'Exp_IronCAA_WT/'
#rootPath=source+'Exp_Incub_4h30-6h30-24h-52h/'
#rootPath=source+'Exp_Incubation_14h-19h-43h/'
#rootPath= source +'Exp_Incubation_merge/'
#rootPath=source+'Exp_DilutionWTstock2/'
rootPath=source+'Exp_Dilution_5_50_500_3500/'
#rootPath=source+'Exp_IronCAA_PvdS/'
#rootPath=source+'Exp_Dilution_PvdS/'
#rootPath=source+'Exp_BipyCAA_WT/'
#rootPath=source+'Exp_BipyCAA_PvdS/'
#rootPath= source +'Exp_IronPvd/'
#rootPath= source +'Exp_GreenRed/'
#rootPath= source +'Exp_Incubation_4h30_6h30_24h_43h13_1bactPerDrop/'
#rootPath='/Users/admin/Documents/Experiences/Milli/Exp_Dilution_M9glycerol/'


rootPathList=[rootPath]
#rootPathList=[source+'Exp_Dilution_PvdS/',\
#              source+'Exp_IronCAA_WT/', source+'Exp_Incub_4h30-6h30-24h-52h/', \
#              source+'Exp_Incubation_14h-19h-43h/',\
#              source+'Exp_Dilution_5_50_500_3500/',source+'Exp_DilutionWTstock2/']
              #source+'Exp_IronCAA_PvdS/',\
#              source +'Exp_Incubation_4h30_6h30_24h_43h13_1bactPerDrop/',source+'Exp_BipyCAA_WT/',source+'Exp_BipyCAA_PvdS/']

print(rootPathList)


# In[59]:


#label=['SBW25-WT_CAA_1000','SBW25-WT_CAA_7000']#'SBW25-WT_CAA_10', 'SBW25-WT_CAA_100', 
for rootPath in rootPathList:
    folder=sorted([join(join(rootPath,o),'analysis/') for o in listdir(rootPath) if (o[:3]=='201' or o[:3]=='EXP') and os.path.isdir(os.path.join(rootPath, o))])
    
    #channelList =['PVD','GFP','RFP','PVDoverGFP', 'SpeedSize']
    channelList =['GFP']
    incCarte=5e-3 #incertitude de la mesure par la carte hardware
    timeThresh=30 #heures max d'analyse
    threshOutlayer=25 #heures pour mesure de contamination. Pour desactiver mettre une val negative
    display = False
    startTime=0 #time to start analysing the curves to get ride of the weird shit at the beginning
    deletionData=False #delete data in resultIndiv
    
    for channel in channelList:
        for path in [folder[2]]:
            print(path)
            [dropMap, df, label]=loadData(path)
            checkfitData(df, dropMap, label, channel,path,startTime,timeThresh,threshOutlayer, incCarte)
            fitDataIndiv(df, dropMap, label, channel,path,startTime, timeThresh, threshOutlayer, incCarte, display, deletionData)
            
            
####not used anymore
#poolDataInterpolate(dropMap, df, label, channel,path, incCarte, nbReps,startTime, timeThresh, threshOutlayer, display)
#checkDataInterpolated(path, label,channel)
#fitDataLM(path, label,channel) #Not used after lukas discussion


# In[62]:


#plot le diagrame a moustache du lag, gRate, yield en allant chercher les bons label directement et en les ordonnant
#sauve les plot et les grandeurs mesurées dans plotIndiv
#save of the single charte npy file

for rootPath in rootPathList:
    folder=[join(join(rootPath,o),'analysis/') for o in listdir(rootPath) if (o[:3]=='201' or o[:3]=='EXP') and os.path.isdir(os.path.join(rootPath, o)) ]


    listParam=['stdgRate','gRate', 'lag', 'yld']
    #listParam=['yld']
    channel='GFP'  #GFP RFP PVD

    for kindOfData in listParam:

        for path in folder:
            print(path)

            a=path.split('/')
            date=a[-3]
            data=[]
            cleanedData=[]
            listData=[]
            listLabel=[]

            [nn, nnn, label]=loadData(path)
            [stdgRate, gRate, lag, stdlag, yld, stdyield]=getDataHistogram(label,path,channel)
            for l in label:

                #defini le yrange
                if kindOfData == 'gRate':
                    data=gRate[l]
                    yMin=0
                    yMax=2
                    #yMax=1 #GFP
                elif kindOfData =='lag':
                    data=lag[l]
                    yMin=0
                    yMax=30
                elif kindOfData == 'yld':
                    data=yld[l]
                    yMin=0
                    yMax=8
                elif kindOfData == 'stdgRate':
                    data=stdgRate[l]
                    yMin=0
                    yMax=.5
                else :
                    print('do not get kind of data')

                cleanedData = [x for x in data if str(x) != 'nan']
                listData.append(cleanedData)
                listLabel.append(l+' N='+str(len(cleanedData)))


            num = []
            idx = []
            j = 0

            print(label)
            for a in label:
                val=getValueInLabel(a,path)

                num.append(val)
                idx.append(j)
                j = j+1

            #ordonne les label selon une concentration croissant
            zipped = zip(num,idx)
            zz = sorted(zipped,key = lambda t: t[0])

            [nn, idx] = [list(t) for t in zip(*zz)]
            listLabelOrdered = [listLabel[i] for i in idx]
            listDataOrdered = [listData[i] for i in idx]

            plt.boxplot(listDataOrdered)

            plt.title(date+'_'+channel)
            plt.grid()
            plt.ylabel(kindOfData)
            plt.ylim([yMin, yMax])
            
            if rootPath==source+'Exp_Incubation_14h-19h-43h/' and rootPath==source+'Exp_Incub_4h30-6h30-24h-52h/':
                plt.xticks(range(0,60,50),rotation=90)
                plt.xlabel('Time of preculture (h)') 
            else:
                plt.xticks(range(1,len(listLabelOrdered)+1), listLabelOrdered)
                plt.xticks(rotation=90)

            pathPlot = path+'plotIndiv/'
            if not os.path.exists(pathPlot):
                os.makedirs(pathPlot)

            plt.savefig(pathPlot+'plot_'+ kindOfData +'_'+date+'_'+channel+'.pdf',bbox_inches='tight', format='pdf')
            plt.show()


            dataDict=dict(zip(listLabelOrdered,listDataOrdered))

            with open(pathPlot+'result_'+ kindOfData +'_'+ channel+'_'+date+'.csv', "w") as f:
                w = csv.writer(f)
                w.writerow(dataDict.keys())
                w.writerow(dataDict.values())

            # Save dict 
            np.save(pathPlot+ kindOfData+'_'+ channel+'.npy', dataDict) 


# =============================================================================
# # In[141]:
# 
# 
# #plot in one chart the parameters measured as a function of the varying (conc, age, bipy, iron) param for the different experiment. The points correspond to median  5% and 95% bar.
# #extraction from the single charte npy file saved in previous cell
# for rootPath in rootPathList:
#     folder=[join(join(rootPath,o),'analysis/') for o in listdir(rootPath) if (o[:3]=='201' or o[:3]=='EXP') and os.path.isdir(os.path.join(rootPath, o)) ]
# 
#     listParam=['stdgRate','gRate', 'lag', 'yld']
#     channel='GFP' #GFP RFP PVD
#     for kindOfData in listParam:
# 
#         plt.figure()
#         c = matplotlib.cm.rainbow(np.linspace(0, 1, len(folder)))
#         plt.rcParams.update({'font.size': 14})
#         idxColor=0
# 
#    
#         #for bipy
#         if kindOfData=='yld':
#             setZero=0
#         elif  kindOfData=='gRate':
#             setZero=0
#         elif kindOfData=='lag':
#             setZero=0
#         else:
#             setZero=0
#         yMax=0
#         #concMax=3000
#         #plt.plot([0,concMax], [setZero,setZero], '--') #plot line for non growing condition
#         
#         #for iron
#         concMax=500
#         
#         # Load
#         for path in folder:
#             print(path)
#             read_dictionary = np.load(path+'plotIndiv/'+kindOfData+'_'+ channel+'.npy',allow_pickle = True).item()
# 
#             conc=[]
#             med=[]
#             perc75=[]
#             perc25=[]
#             stdv=[]
#             
#             
#             for key in read_dictionary.keys():
#                 tmp = np.array(read_dictionary[key])
#                 val=getValueInLabel(key,path)
# 
#                 
#                 if len(tmp)>5:
# 
#                     conc.append(val)
#                     tmpMed=np.percentile(tmp, 50)
#                     med.append(tmpMed) #median
#                     perc25.append(tmpMed-np.percentile(tmp, 25))
#                     perc75.append(np.percentile(tmp, 75)-tmpMed) 
#                     stdv.append(np.std(tmp))
#                 else:
#                     #if did not grow in more than 5 drops over 80 set to zero
#                     conc.append(val)
#                     med.append(setZero) #median
#                     perc75.append(0) #95% bar
#                     perc25.append(0) #5% bar
#                     stdv.append(0)
# 
#             if yMax<np.max(np.array(med)+np.array(perc75)):
#                 yMax=np.max(np.array(med)+np.array(perc75))
#             #plot for every exp with different color
# 
#             plt.errorbar(conc, med, yerr=[perc25,perc75], fmt='x', color=c[idxColor], label=str(idxColor),uplims=True, lolims=True, markersize=10)
#             idxColor=idxColor+1
#             
#         date=rootPath.split('/')
# 
#         plt.title(date[-2]+' '+kindOfData+' '+channel)
#         #plt.ylim(setZero, yMax*1.1)
# 
#         plt.grid()
#         #plt.legend()
#         
#         axis_font = {'fontname':'Courrier', 'size':'14'}
#         if kindOfData=='yld':
#             plt.ylabel('yield (au)',**axis_font)
#             plt.ylim([-5,0])
#         elif  kindOfData=='gRate':
#             plt.ylabel('growth rate (1/h)',**axis_font)
#             plt.ylim([0.7,.9])
#         elif kindOfData=='lag':
#             plt.ylabel('lag (h)',**axis_font)
#             plt.ylim([0,30])
#         elif kindOfData=='stdgRate':
#             plt.ylabel('error gRate (1/h)',**axis_font)
#             plt.ylim([0,.1])
# 
#         
#         #manip dilution
#         #plt.xscale('log')
#         #plt.xlabel('inoculum (cell/drop)',**axis_font)
#         #manip iron
#         #plt.xticks(rotation=90)
#         #plt.xlabel('Pre-incubation time (h)')
#         #manip iron
#         #plt.xticks(range(0,concMax,100),rotation=90)
#         #plt.xlabel('Fe2(SO4)3 (µM)')
#         #manip bipy
#         #plt.xticks(range(0,3000,200),rotation=90)
#         #plt.xlabel('Bipyridyl (µM)')
#         if rootPath==source+'Exp_GreenRed/':
#                 plt.xticks(rotation=90)
#             
#         if rootPath==source+'Exp_Incubation_14h-19h-43h/' or rootPath==source+'Exp_Incub_4h30-6h30-24h-52h/' or rootPath==source+'Exp_Incubation_merge/':
#             plt.xticks(range(0,60,10))
#             plt.xlabel('Age of preculture (h)') 
#             
#         plt.savefig(rootPath+kindOfData+'_'+ channel+'.pdf',bbox_inches='tight', format='pdf')
#         plt.savefig(rootPath+kindOfData+'_'+ channel+'.jpeg',bbox_inches='tight', format='jpeg')
#         plt.show()
# 
# =============================================================================

# =============================================================================
# # In[75]:
# 
# 
# #plot correlation between lag and gRate 
#  
# clr = matplotlib.cm.rainbow(np.linspace(0, 1, 12))
# for rootPath in rootPathList:
#     folder=[join(join(rootPath,o),'analysis/') for o in listdir(rootPath) if (o[:3]=='201' or o[:3]=='EXP') and os.path.isdir(os.path.join(rootPath, o)) ]
# 
# 
#     channel='GFP'  #GFP RFP PVD
# 
# 
#     
#     for path in folder:
#         print(path)
# 
#         a=path.split('/')
#         date=a[-3]
# 
#         idxClr=0
# 
#         [nn, nnn, label]=loadData(path)
#         [stdgRate, gRate, lag, stdlag, yld, stdyield]=getDataHistogram(label,path,channel)
# 
# 
#         plt.figure()
#         for l in label:
# 
# 
#             plt.scatter(gRate[l],lag[l], color=clr[idxClr])
#             idxClr=1+idxClr
#             plt.ylabel('lag')
#             plt.xlabel('gRate')
#             
#         plt.grid()
#         plt.xlim([0,1.2])
#         plt.ylim([3,30])
#         plt.savefig(path+'correlLagGrate_'+kindOfData+'_'+ channel+'.pdf',bbox_inches='tight', format='pdf')
#         plt.show()
# 
# =============================================================================

# In[]
#plot CUMULATIVE distribution of gRate, lag, yield in a facet for each conditions
#extraction from the single charte npy file saved in previous cell
for rootPath in rootPathList:
    folder=[join(join(rootPath,o),'analysis/') for o in listdir(rootPath) if (o[:3]=='201' or o[:3]=='EXP') and os.path.isdir(os.path.join(rootPath, o))]


    listParam=['stdgRate','gRate', 'lag', 'yld']
    

    channel='GFP' #GFP RFP PVD
    for kindOfData in listParam:
        
        fig, ax = plt.subplots(1, 1)

        for m_path in folder:
            print(m_path)

            read_dictionary = np.load(m_path+'plotIndiv/'+kindOfData+'_'+ channel+'.npy',allow_pickle = True).item()
    
            df=pddf(dict([ (k,pd.Series(v)) for k,v in read_dictionary.items() ])) #cree un dataframe a partir du dictionnaire
            
            for key in df.keys(): 
                k=key.split(' ')
                label=getValueInLabel(k[0],m_path)
                print([key, label])
                df = df.rename(columns={key: label})           #change le nom des columns pour enlever le N=... 
                
                
            for key in df.keys():    
                
                data=df[key][np.isfinite(df[key])] #remove the nan values
                # evaluate the histogram
                values, base = np.histogram(data, bins='fd')
                #evaluate the cumulative
                cumulative = np.cumsum(values)

                

                if key==5 or key==4.5 or key==4 or key==14 or key==1:
                    if '2018-02-21_dilution_CAA_WT-stock2' in m_path:
                        ax.plot(base[:-1], cumulative/len(data), '-*k', markerSize=10)
                    else:
                        ax.plot(base[:-1], cumulative/len(data), '-k')
                elif key==50 or key==25 or key==6.5 or key==6 or key==10:
                    ax.plot(base[:-1], cumulative/len(data), '-b')
                elif key==500 or key==100 or key==24  or key==19 or key==100:
                    ax.plot(base[:-1], cumulative/len(data), '-r')
                elif key==3500 or key==1000 or key==52 or key==43 or key==1000:
                    ax.plot(base[:-1], cumulative/len(data), '-g')                                    
                else :
                    print('not printed:'+str(key))
                
            
        lgd=my_legend(rootPath)
            
        if kindOfData=='gRate':
            xl=0.3
            xll=1.5
            plt.xlim([xl,xll])
            plt.xticks(np.arange(xl,xll,.05),rotation=60)
            plt.text(xll*0.8, 0.5,lgd) 
            
        if kindOfData=='yld':
            xl=5
            xll=8
            plt.xlim([xl,xll])
            plt.xticks(np.arange(xl,xll,.5))
            plt.text(xll*0.9, 0.5,lgd) 

        if kindOfData=='lag':
            xl=0
            xll=20
            plt.xlim([xl,xll])
            plt.xticks(np.arange(xl,xll,1),rotation=60)
            plt.text(xll*0.8, 0.1,lgd) 



        plt.yticks(np.arange(0,1.1,.1))
        plt.ylim([0,1.1])            

        plt.xlabel(kindOfData)
        plt.ylabel('cumulative frequency of '+kindOfData)
        plt.box   
        plt.grid(which='major', axis='both')
        plt.savefig(rootPath+kindOfData+'_cumul'+ channel+'.pdf',bbox_inches='tight', format='pdf')
        plt.savefig(rootPath+kindOfData+'_cumul'+ channel+'.jpeg',bbox_inches='tight', format='jpeg')
        plt.show()
  
# In[]
#plot HISTO distribution of gRate, lag, yield in a facet for each conditions
#extraction from the single charte npy file saved in previous cell
for rootPath in rootPathList:
    folder=[join(join(rootPath,o),'analysis/') for o in listdir(rootPath) if (o[:3]=='201' or o[:3]=='EXP') and os.path.isdir(os.path.join(rootPath, o))]


    listParam=['gRate']#['stdgRate','gRate', 'lag', 'yld']
    

    channel='GFP' #GFP RFP PVD
    for kindOfData in listParam:
        
        fig, ax = plt.subplots(1, 1)

        for m_path in folder:
            print(m_path)

            read_dictionary = np.load(m_path+'plotIndiv/'+kindOfData+'_'+ channel+'.npy',allow_pickle = True).item()
    
            df=pddf(dict([ (k,pd.Series(v)) for k,v in read_dictionary.items() ])) #cree un dataframe a partir du dictionnaire
            
            for key in df.keys():
                k=key.split(' ')
                label=getValueInLabel(k[0],m_path)
                print([key, label])
                df = df.rename(columns={key: label})           #change le nom des columns pour enlever le N=... 
                
                
            for key in df.keys():    
                
                data=df[key][np.isfinite(df[key])] #remove the nan values
                # evaluate the histogram
                values, base = np.histogram(data, bins='fd')
                #evaluate the cumulative


               
                plot_distribution(ax,key,base,values/len(values),m_path)
# =============================================================================
# 
#                 if key==5 or key==4.5 or key==4:
#                     if '2018-02-21_dilution_CAA_WT-stock2' in m_path:
#                         ax.plot(base[:-1], values/len(values), '-*k', markerSize=10)
#                     else:
#                         ax.plot(base[:-1], values/len(values), '-ok')
#                 elif key==50 or key==25 or key==6.5 or key==6:
#                     ax.plot(base[:-1], values/len(values), '-ob')
#                 elif key==500 or key==100 or key==24:
#                     ax.plot(base[:-1], values/len(values), '-or')
#                 elif key==3500 or key==1000 or key==52 or key==43:
#                     ax.plot(base[:-1], values/len(values), '-og')                                    
#                 else :
#                     print('not printed:'+str(key))
# =============================================================================
        
        lgd=my_legend(rootPath)
            
        if kindOfData=='gRate':
            plt.xlim([0.6,0.9])
            plt.xticks(np.arange(0.6,1,.02),rotation=60)
            plt.text(0.61, 2,lgd) 
            
        if kindOfData=='yld':
            plt.xlim([6,9])
            plt.xticks(np.arange(6,9,.5))
            plt.text(6.5, 2,lgd) 

        if kindOfData=='lag':
            plt.xlim([0,20])
            plt.xticks(np.arange(0,20,1),rotation=60)
            plt.text(13, 2,lgd) 


        plt.yticks(np.arange(0,4))
        plt.ylim([0,4])            

        plt.xlabel(kindOfData)
        plt.ylabel('cumulative frequency of '+kindOfData)
        plt.box   
        plt.grid(which='major', axis='both')
        plt.savefig(rootPath+kindOfData+'_histo'+ channel+'.pdf',bbox_inches='tight', format='pdf')
        plt.savefig(rootPath+kindOfData+'_histo'+ channel+'.png',bbox_inches='tight', format='png')
        plt.show()
        








# In[]

        

ax=df[3500.0].plot.hist(color='b',edgecolor='black', linewidth=1.2)
ax.set_xlabel('growth rate')
ax.set_ylabel('frequency')
ax.set_xlim([0.78, 0.92])
plt.savefig(rootPath+kindOfData+'_histo3500'+ channel+'.png',bbox_inches='tight', format='png')



# In[]
        
#plot distribution of gRate, lag, yield in a facet for each conditions
#extraction from the single charte npy file saved in previous cell
for rootPath in rootPathList:
    folder=[join(join(rootPath,o),'analysis/') for o in listdir(rootPath) if (o[:3]=='201' or o[:3]=='EXP') and os.path.isdir(os.path.join(rootPath, o))]


    listParam=['stdgRate','gRate', 'lag', 'yld']
    

    channel='GFP' #GFP RFP PVD
    for kindOfData in listParam:
        
        df_concat=pddf()
        rep=0
        #eqRepDate=pddf()
        
    

        for m_path in folder:
            print(m_path)
            print(rep)
            read_dictionary = np.load(m_path+'plotIndiv/'+kindOfData+'_'+ channel+'.npy',allow_pickle = True).item()
    
            df=pddf(dict([ (k,pd.Series(v)) for k,v in read_dictionary.items() ])) #cree un dataframe a partir du dictionnaire
            
            for key in df.keys():
                k=key.split(' ')
                label=getValueInLabel(k[0],m_path)
                print([key, label])
                df = df.rename(columns={key: label})           #change le nom des columns pour enlever le N=... 


            df=df.join(pddf(data=np.full(len(df.index), str(m_path.split('/')[-3].split('_')[0])),columns=['rep']))   #ajoute une column numero pour la manip (jour)      
            
            #eqRepDate=eqRepDate.append(pddf({'file':m_path.split('/')[-3],'rep':pd.Series(rep)}), ignore_index=True)



            rep+=1 
            
            if rep==0:
                df_concat=df
            else:
                df_concat=df_concat.append(df)     #stock toute les experiences dans un même dataframe


        #eqRepDate.to_csv(rootPath+'Distri_eqRepToDate.csv')
        
        if df_concat.size>0:
            #supprime des columns mauvaise
            if rootPath=='/Users/admin/Documents/Experiences/Milli/Exp_Incubation_4h30_6h30_24h_43h13_1bactPerDrop/':
                df_concat=df_concat.drop([0.0, 6.0], axis=1)     
        
        
            if rootPath=='/Users/admin/Documents/Experiences/Milli/Exp_Dilution_5_50_500_3500/':
                df_concat=df_concat.drop([0.5, 350], axis=1) 
            
            df_melt=pddf()
            df_melt=pd.melt(df_concat, id_vars="rep") #transform le dataframe avec les replicats en colomuns en un dataframe avec une trois columns: nom de la condition(entête des columns), numero de manip, valeur
            df_melt=df_melt[np.isfinite(df_melt['value'])] #remove the nan values
        
        
         
            if kindOfData=='gRate':
                xlow = 0
                xhight = 1.3
                xlabel = 'growth rate (1/h)'
            
            if kindOfData== 'stdgRate':
                xlow = 0
                xhight = .1
                xlabel = 'error growth rate (1/h)'
        
                
            elif kindOfData=='yld':
                xhight = -1
                xlow = -3
                xlabel = 'yield (au)'
                
            elif kindOfData=='lag':
                xlabel = 'lag (h)'
                xlow = 3
                xhight = 30
                
            else :
                xlow = None
                xhight = None
                xlabel = kindOfData
        
        
            if 'p' in locals():
                del p
                
            if len(df_melt['value']) > 0:
                bwth=2 * scipy.stats.iqr(df_melt['value']) / (len(df_melt['value'])**(1/3))
            else:
                bwth=None
            
        #geom_density(bw=1/100) \
          #  p=ggplot(aes(x='value', fill='rep', colour='rep'), data=df_melt) \
          #        +  geom_histogram(binwidth=bwth, alpha=.5) \
          #        +  facet_wrap('variable') \
          #        +  xlim(low=xlow, high=xhight) \
          #        +  ylab('occurence') \
          #        +  xlab(xlabel) \
         
                  
        
        
                 
                
            m_name = rootPath+'distri_' + kindOfData + '_'+ channel +'.pdf'
            #p.save(m_name)
            
            




# In[405]:


#plot the difference of lag of PVD and GFP.

kindOfData='lag'
setZero=0
dfLag=pddf()
dfLagStd=pddf()
clm=0
flag=True

for rootPath in rootPathList:
    folder=[join(join(rootPath,o),'analysis/') for o in listdir(rootPath) if o[:3]=='201' or o[:3]=='EXP']

    for path in folder:
        print(path)

        concGFP=[]
        medGFP=[]
        medPVD=[]
        concPVD=[]
        stdvGFP=[]
        stdvPVD=[]
        
        channel='GFP'
        read_dictionary = np.load(path+'plotIndiv/'+kindOfData+'_'+ channel+'.npy').item()

        for key in read_dictionary.keys():
            tmp = np.array(read_dictionary[key])
            print(key)
            val=getValueInLabel(key,path)


            if len(tmp)>5:

                concGFP.append(val)
                tmpMed=np.mean(tmp)
                medGFP.append(tmpMed) #median
                stdvGFP.append(np.std(tmp)/len(tmp))

            else:
                #if did not grow in more than 5 drops over 80 set to zero
                concGFP.append(val)
                medGFP.append(setZero) #median
                stdvGFP.append(0)

        channel='PVD'
        read_dictionary = np.load(path+'plotIndiv/'+kindOfData+'_'+ channel+'.npy').item()


        for key in read_dictionary.keys():
            tmp = np.array(read_dictionary[key])
            val=getValueInLabel(key,path)


            if len(tmp)>5:

                concPVD.append(val)
                tmpMed=np.mean(tmp)
                medPVD.append(tmpMed) #median
                stdvPVD.append(np.std(tmp)/np.sqrt(len(tmp)))

            else:
                #if did not grow in more than 5 drops over 80 set to zero
                concPVD.append(val)
                medPVD.append(np.nan) #median
                stdvPVD.append(np.nan)

        
        
        if flag:
            flag=False
        else:
            clm=clm+1
            sLag=pd.Series(np.array(medPVD)-np.array(medGFP),  index=np.array(concPVD))
            dfLag[clm]=sLag
            sLag=pd.Series(np.array(stdvPVD)+np.array(stdvGFP),  index=np.array(concPVD))
            dfLagStd[clm]=sLag

print(dfLagStd)
m=dfLag.mean(axis = 1)
s=dfLagStd.mean(axis = 1)

nbsample=len(dfLagStd.columns.values)


plt.figure()
plt.errorbar(m.index.values,m.values,s.values,fmt='x',markersize=20,c='k')
plt.grid()
plt.ylabel('lagPVD-lagGFP (h)')
plt.xlabel('bacterial concentration (cell/drop)')
plt.ylim([0,6])
plt.xscale('log')

plt.savefig(rootPath+'diffLag'+ '.pdf',bbox_inches='tight', format='pdf')
plt.show()


# In[58]:


#get the diversity of pvd expression per bact within a drop as a function of the controle parameters of the exp

folder=sorted([join(join(rootPath,o),'analysis/') for o in listdir(rootPath) if o[:3]=='201' or o[:3]=='EXP'])
channel='PVDoverGFP'
timeThresh=35 #heures max d'analyse
threshOutlayer=20#heures pour mesure de contamination
startTime=6
seuilProd=0.5
for path in folder:
    print(path)
    [dropMap, df, label]=loadData(path)

    dfDiv=pddf()
    for k,dataFile in enumerate(label):
        sr = pd.Series()
        tmp= pddf()        
        for j,content in enumerate(dropMap[:,1]):
            
            if content==dataFile:    
                y=[]
                x=[] 
                if dataFile=='CAA-WT-100':
                    startTime=10
                elif dataFile=='CAA-WT-25':
                    startTime=15
                else:
                    startTime=6
                    
                    
                [x,y,flag] = getfitData(df, j, channel,startTime,timeThresh,threshOutlayer) 
                if flag and len(x)>0:
                    idxProd=bisect.bisect(y,seuilProd)
                    sr=sr.append(pd.Series(x[idxProd]),ignore_index=True)
        tmp[dataFile]=sr
        dfDiv = pd.concat([dfDiv,tmp], ignore_index=True)

    

    print(dfDiv.std(axis = 0))
      


# In[730]:


#Exp_Dilution_5_50_500_3500
conc0=5    
val7=1.166799
conc1=np.array([3500,5,50,500])
val5=np.array([0.525821,1.418746, 0.645083, 0.570822])
val6=np.array([0.430209,1.486276,0.568323,0.388596])

#Exp_DilutionWTstock2
conc2=np.array([100,1000,25])
val1=np.array([1.008576, 0.530170, -10])
val2=np.array([0.566269, 0.466830, 1.089643])
val3=np.array([0.643435, 0.561668, 1.051770])
val4=np.array([0.641859, 0.646762, 1.700297])


rootPath==source+'Exp_Dilution_5_50_500_3500/'
if False:
    plt.figure()
    plt.scatter(conc0,val7)
    plt.scatter(conc1,val5)
    plt.scatter(conc1,val6)
    plt.xscale('log')
    plt.grid('on')
    plt.ylim([0,2])
    plt.xlabel('bacteria inoculum (cell/drop)')
    plt.ylabel('std production of pvd per cell')
    plt.savefig(rootPath+'diversityPVD'+ '.pdf',bbox_inches='tight', format='pdf')
    plt.show()

rootPath==source+'Exp_DilutionWTstock2/'


if False:
    plt.figure()
    plt.scatter(conc2,val1)
    plt.scatter(conc2,val2)
    plt.scatter(conc2,val3)
    plt.scatter(conc2,val4)
    plt.xscale('log')
    plt.grid('on')
    plt.ylim([0,2])
    plt.xlabel('bacteria inoculum (cell/drop)')
    plt.ylabel('std production of pvd per cell')
    plt.savefig(rootPath+'diversityPVD'+ '.pdf',bbox_inches='tight', format='pdf')
    plt.show()
   
#All
plt.figure()
plt.scatter(conc0,val7)
plt.scatter(conc1,val5)
plt.scatter(conc1,val6)
plt.scatter(conc2,val1)
plt.scatter(conc2,val2)
plt.scatter(conc2,val3)
plt.scatter(conc2,val4)

plt.xscale('log')
plt.grid('on')
plt.ylim([0,2])
plt.xlim([-1,4000])
plt.xlabel('bacteria inoculum (cell/drop)')
plt.ylabel("std production \nof pvd per cell (h)")
plt.savefig(rootPath+'diversityPVD_all'+ '.pdf',bbox_inches='tight', format='pdf')
plt.title('all')
plt.show()



# In[60]:


#plot the probability of growing as a function of bipy concentration
channel='GFP' #GFP RFP
for rootPath in rootPathList:
    folder=[join(join(rootPath,o),'analysis/') for o in listdir(rootPath) if o[:3]=='201' or o[:3]=='EXP']

    plt.figure()
    clr = matplotlib.cm.rainbow(np.linspace(0, 1, len(folder)))
    mkr=['x','+','<','>','.','x','+','<','>','.','x','+','<','>','.']
    idxColor=0

    for path in folder:
        pathPlot = path+'plotIndiv/'
        
        [dropMap, nnn, label]=loadData(path)
        [stdgRate, gRate, lag, stdlag, yld, stdyield]=getDataHistogram(label,path,channel)
         

        dfProba=pddf()
        for l in label:
            data=yld[l]
            cleanedData = [x for x in data if str(x) != 'nan']
            conc=getValueInLabel(l,path)
            nbTot=np.sum(dropMap[:,1]==l)
            dfProba[conc] = pd.Series([len(cleanedData), nbTot, len(cleanedData)/nbTot], index=[0, 1, 2])


        #save number of probability drop
        with open(pathPlot+'result_'+ 'proba' + channel+'.csv', "w") as f:
            w = csv.writer(f)
            #w.writerow(dataN.keys())
            #w.writerow(dataN.values())
            w.writerow(dfProba.columns.values)
            w.writerow(dfProba.iloc[2].values)




        #setZero=setZero-setZero/10 #change the position of the non growing drop on the figure to see when they stack on eacho other

        conc = list(map(int, dfProba.columns.values))#convert list of string to list of int
        proba=dfProba.iloc[2].values

        plt.plot(conc, proba,mkr[idxColor], color=clr[idxColor], label=str(idxColor))
        idxColor=idxColor+1
    date=rootPath.split('/')
    print(date[-2])
    plt.title(date[-2]+' proba of growing')
    plt.ylim(-0.5,1.5)
    plt.grid()
    plt.ylabel('# growing/# total')
    

    if rootPath==source+'Exp_Incubation_14h-19h-43h/' or rootPath==source+'Exp_Incub_4h30-6h30-24h-52h/' or rootPath==source+'Exp_Incubation_merge/':
        plt.xticks(range(0,60,10))
        plt.xlabel('Age of preculture (h)') 

    elif rootPath==source+'Exp_IronCAA_WT/':
            plt.xticks(range(0,500,200),rotation=90)
            plt.xscale('log')
            plt.xlabel('Fe2(SO4)3 (µM)')
    else:
        #manip bipy
        plt.xticks(range(0,3000,200),rotation=90)
        plt.xlabel('Bipyridyl (µM)')
        #plt.legend()

    plt.savefig(rootPath+'proba'+channel+'.pdf',bbox_inches='tight', format='pdf')
    plt.show()



# In[54]:
#define fonction


# In[10]:


gRateList=[]
stdgRateList=[]
lagList=[]
stdlagList=[]
yieldList=[]
stdyieldList=[]
labelList=[]

channel='GFP'

#WORKS FOR fitDataLM!
for path in folder:
    #print(path)    
    #path='/Users/maxime/Documents/experiences/milli/expDilution_5_50_500_3500/2018-02-09_dilution_CAA_WT-stock2/analysis/'

    fileList=glob.glob(path+'*'+channel+'resultInterp.csv')
    
    
    for file in fileList:
        with open(file, 'r') as f:
            reader = csv.reader(f)
            your_list = list(reader)
        
        name=os.path.split(file)
        #print(your_list)

        
        stdgRateList.append(float(your_list[findIdx(your_list,'max df std')][1]))
        gRateList.append(float(your_list[findIdx(your_list,'max df')][1]))
        lagList.append(float(your_list[findIdx(your_list,'lag time')][1]))
        stdlagList.append(float(your_list[findIdx(your_list,'lag time std')][1]))
        yieldList.append(float(your_list[findIdx(your_list,'max y')][1]))
        stdyieldList.append(float(your_list[findIdx(your_list,'max y std')][1]))
        labelList.append(name[1][:-16])


# In[11]:


# In[ ]:





# In[12]:


def fitDataLM(folder,label,channel): 
    
    for j,dataFile in enumerate(label):
        print(dataFile)
        file=folder+dataFile+channel+'Interp.csv'
        print(file)
        #data=pd.read_csv(file)
        # load data        
        data= np.loadtxt(file,skiprows=1, delimiter=',',dtype=float)
        #od=data.values[:,:-2]

        try:
            t, od= data[:,0], data[:, 1:]
            #print(od)

            # redefine bounds and run inference
            #hyperparameters 
            #0: amplitude
            #1: flexibility
            #2: measurement error
            #initial    b= {0: [-1,5], 1: [-5,2], 2: [-7,1]}
            b= {0: [-1,5], 1: [-5,0], 2: [-1,0]}
            q= fitderiv(t, od, bd= b, logs= False,esterrs=True)

            # plot results
            fig=plt.figure()
            #plt.subplot(2,1,1)
            q.plotfit('f')
            plt.title(dataFile)
            plt.show()
            fig.savefig(folder+dataFile+'fitInterp'+channel+'.jpg')
            #plt.subplot(2,1,2)
            fig=plt.figure()
            q.plotfit('df')
            plt.title(dataFile)
            plt.show()
            fig.savefig(folder+dataFile+'derivativeInterp'+channel+'.jpg')


            # export results
            statdict=q.printstats()
            with open(folder+dataFile+channel+'resultInterp.csv','w') as csv_file:
                writer = csv.writer(csv_file)
                for key, value in statdict.items():
                    writer.writerow([key, value])
                    
        except Exception as inst:
            print('bug od')


# In[499]:




#orderList=['300uMGFP','500uMGFP','700uMGFP','900uMGFP','1100uMGFP','1300uMGFP','1500uMGFP','1700uMGFP','1900uMGFP','2100uMGFP','2300uMGFP','2500uMGFP']
#orderList=['300uMPVD','500uMPVD','700uMPVD','900uMPVD','1100uMPVD','1300uMPVD','1500uMPVD','1700uMPVD','1900uMPVD','2100uMPVD','2300uMPVD','2500uMPVD']


#orderList=['0uMGFP','10uMGFP','50uMGFP','100uMGFP','200uMGFP','300uMGFP','400uMGFP','500uMGFP','1000uMGFP','2000uMGFP','2500uMGFP','3000uMGFP']
#orderList=['0uMPVD','10uMPVD','50uMPVD','100uMPVD','200uMPVD','300uMPVD','400uMPVD','500uMPVD','1000uMPVD','2000uMPVD','2500uMPVD','3000uMPVD']

#orderList=['0uMGFP','10uMGFP','50uMGFP','100uMGFP','200uMGFP','300uMGFP','400uMGFP','500uMGFP','700uMGFP','900uMGFP','1000uMGFP','1200uMGFP']
#orderList=['0uMPVD','10uMPVD','50uMPVD','100uMPVD','200uMPVD','300uMPVD','400uMPVD','500uMPVD','700uMPVD','900uMPVD','1000uMPVD','1200uMPVD']


orderList=['CAA-WT-0uMGFP', 'CAA-WT-5uMGFP', 'CAA-WT-10uMGFP', 'CAA-WT-20uMGFP', 'CAA-WT-40uMGFP', 'CAA-WT-60uMGFP', 'CAA-WT-80uMGFP', 'CAA-WT-120uMGFP', 'CAA-WT-170uMGFP', 'CAA-WT-200uMGFP', 'CAA-WT-300uMGFP', 'CAA-WT-500uMGFP']


order=[]
order=getOrderLabel(orderList,labelList)
labelList=organizeLabel(order,labelList)
stdgRateList=organizeLabel(order,stdgRateList)
gRateList=organizeLabel(order,gRateList)
#lagList=organizeLabel(order,lagList)
#stdlagList=organizeLabel(order,stdlagList)
#yieldList=organizeLabel(order,yieldList)
#stdyieldList=organizeLabel(order,stdyieldList)

barPlot(labelList,gRateList,stdgRateList,0,.8,'Growth rate','growth rate ($ h^{-1}$)',rootPath+'gRateInterp.png')
barPlot(labelList,lagList,stdlagList,0,30,'Lag','Lag (h)',rootPath+'lagInterp.png')
barPlot(labelList,yieldList,stdyieldList,-10,0,'Yield','Yield',rootPath+'yieldInterp.png')

