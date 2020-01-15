#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 20 14:35:52 2019

@author: bacterie
"""
import seaborn as sn
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

        

# In[]
path="/Users/bacterie/Documents/Experiences/Milli/EXP_testSort/"
#folder="EXP20191218_1555/"
folder= 'EXP20200109_1646/'

my_df=pd.read_csv(path+folder+"SortMap.dat", sep=' ')
w=my_df['Well']
line = pd.Series([l[0] for l in w ])
col = pd.Series([l[1:] for l in w ])
col = pd.to_numeric(col)


my_df['Line'] = line
my_df['Col'] = col


my_df = my_df.sort_values(by = ['Line','Col'])


#####create a 2D dataframe with the same configuation than the plate filled by a id_cumber of the drops in each well
my_label=pd.DataFrame()
#my_label['label'] = my_df['Well'] + " " + my_df['Waste/Sample'] + " " + my_df['StartDrop'].map(str) +"-"+my_df['EndDrop'].map(str)
my_label['label'] = my_df['Waste/Sample'] + " " + my_df['StartDrop'].map(str) +"-"+my_df['EndDrop'].map(str)
my_label['Col'] = my_df['Col']
my_label['Line'] = my_df['Line']

my_line = my_label.groupby('Line')['label'].apply(list)
my_line['A'].append('A12 Null')
my_line = pd.DataFrame(my_line)

my_plate = pd.DataFrame(my_line['label'].values.tolist(), columns=[i+1 for i in range(12)]) #converte the list in the column 'label' as columns in a dataframe
my_plate.index=my_line.index




#####create a 2D dataframe with the same configuation than the plate filled by a label 1 if the well containe a sample and 0 if the well contain a waste
my_sample=pd.DataFrame()
my_sample['label'] = my_df['Waste/Sample']
my_sample['Col'] = my_df['Col']
my_sample['Line'] = my_df['Line']

my_sample_line = my_sample.groupby('Line')['label'].apply(list)
my_sample_line['A'].append('W')
my_sample_line = pd.DataFrame(my_sample_line)

my_sample_plate = pd.DataFrame(my_sample_line['label'].values.tolist(), columns=[i+1 for i in range(12)]) #converte the list in the column 'label' as columns in a dataframe
my_sample_plate.index=my_sample_line.index
my_sample_plate = my_sample_plate.replace('S', 1)
my_sample_plate = my_sample_plate.replace('W', 0)



#####create a 2D dataframe with the same configuation than the plate filled by the number of drops in each well
my_count=pd.DataFrame()
my_count['label'] = my_df['EndDrop'] - my_df['StartDrop'] + 1
my_count['Col'] = my_df['Col']
my_count['Line'] = my_df['Line']

my_count_line = my_count.groupby('Line')['label'].apply(list)
my_count_line['A'].append(0)
my_count_line = pd.DataFrame(my_count_line)

my_count_plate = pd.DataFrame(my_count_line['label'].values.tolist(), columns=[i+1 for i in range(12)]) #converte the list in the column 'label' as columns in a dataframe
my_count_plate.index=my_count_line.index


my_count_plate = pd.DataFrame(my_sample_plate*my_count_plate, columns=my_count_plate.columns, index=my_count_plate.index) #set to zero the wells of waste. The other wells are labeled by the number of drop



#ax = plt.subplot(111)
#ax = plt.subplot(111, frame_on=False) # no visible frame
#ax.xaxis.set_visible(False)  # hide the x axis
#ax.yaxis.set_visible(False)  # hide the y axis

#pd.plotting.table(ax, data = my_plate)  # where df is your data frame


#plt.savefig('mytable.png')

heatmap_binary(my_count_plate,my_plate)


# In[]:

##
import matplotlib.colors as mcolors
from matplotlib.colors import LinearSegmentedColormap

def heatmap_binary(df, df_label ,edgecolors='w',log=False, path='/Users/bacterie/Desktop/'):    
    width = len(df.columns)/7*10
    height = len(df.index)/7*10
    
    fig, ax = plt.subplots(figsize=(20,10))#(figsize=(width,height))
    

    colors = ['black', 'blue','green']

    cmap = LinearSegmentedColormap.from_list('name', colors)
    norm = plt.Normalize(0, 5)
    
    heatmap = ax.pcolor(df ,
                        edgecolors=edgecolors,  # put white lines between squares in heatmap
                        cmap=cmap,
                        norm=norm,
                        linewidths=1)
    data = df_label.values
    for y in range(data.shape[0]):
        for x in range(data.shape[1]):
            plt.text(x + 0.5 , y + 0.5, '%s' % data[y, x], #data[y,x] +0.05 , data[y,x] + 0.05
                 horizontalalignment='center',
                 verticalalignment='center',
                 color='w',
                 fontsize=10)
    
    
    ax.autoscale(tight=True)  # get rid of whitespace in margins of heatmap
    ax.set_aspect('equal')  # ensure heatmap cells are square
    ax.xaxis.set_ticks_position('top')  # put column labels at the top
    ax.tick_params(bottom='off', top='off', left='off', right='off')  # turn off ticks
    
    ax.set_yticks(np.arange(len(df.index)) + 0.5)
    ax.set_yticklabels(df.index, size=20)
    ax.set_xticks(np.arange(len(df.columns)) + 0.5)
    ax.set_xticklabels(df.columns, rotation=90, size= 15)
    
    # ugliness from http://matplotlib.org/users/tight_layout_guide.html
    #from mpl_toolkits.axes_grid1 import make_axes_locatable
    #divider = make_axes_locatable(ax)
    #cax = divider.append_axes("right", "3%", pad="1%")
    #fig.colorbar(heatmap, cax=cax)
    plt.gca().invert_yaxis()
    fig.savefig(path+'sortedPlate.jpeg', format='jpeg', dpi=800)
                 

#df1 = pd.DataFrame(np.random.choice([0, 0.75], size=(4,5)), columns=list('ABCDE'), index=list('WXYZ'))
#heatmap_binary(df1)
    
# In[]:
for i in range(49,60,2):
    print(str(i)+"\n")
for i in range(271,290,2):
    print(str(i)+"\n")
for i in range(501,520,2):
    print(str(i)+"\n")
for i in range(731,750,2):
    print(str(i)+"\n")

                 
                 
                