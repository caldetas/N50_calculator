#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 29 14:34:57 2018

@author: hannes
"""


class data:
    def __init__(self):
            self.df = {}
            self.settings = {}
    def add_df(self, genome, db, tipe, df):
        self.df[tipe][db][genome] = df
    def add_setting(self, name, settings):
        self.settings[name] = settings

#init class data
data = data()


import sys
import getopt
import time


#instructions
form='N50.py\n\t--minc *minimum contig length, opt*\n\t <flat>'

data.add_setting('time1', time.time())

def main(argv):

    minc = 0
    
    try:
        opts, args = getopt.getopt(argv,"h",["minc="])
    except getopt.GetoptError:
        print ('{}'.format(form))
        sys.exit()
    for opt, arg in opts:        
        if opt == '-h' or opt == '-help' or opt == '--help':
            print ('{}'.format(form))
            sys.exit()
        elif opt == '--minc':
            minc = int(arg)

    print()
#    print('argv =', argv)
    
    print()
    
      
    import numpy as np
    import re
    import pandas as pd

    l1=[]
    t=0
    GC_count = 0

    # open the fasta file
    ass=open('{}'.format(argv[-1]))
    for line in ass:
        line = line.strip('\n')

        #check if .fasta header
        match= re.search("^(>)", line)

        #first contig  
        if match and t == 0:
            l1.append([line])           
            temp=''
            t=1
        #if not first contig
        elif match and t==1:

            #length of sequence
            l1[-1].append(len(temp))
            l1.append([line])
            if len(temp) >= minc:
                GC_count += temp.count('G') + temp.count('C')
            temp=''
        else:
            temp += line

    #last contig
    l1[-1].append(len(temp))
    if len(temp) >= minc:
        GC_count += temp.count('G') + temp.count('C')
    temp=''
    t=0

    ass.close()
   
    data.add_setting('time2', time.time())
    
    #create pandas dataframe of fasta file
    df = pd.DataFrame(l1)
    df = df.iloc[:,:2]
    df = df.sort_values(by=[1])
    df = df.reset_index(drop=True)
    df = df.loc[df[1] >= minc]
    df = df.rename(columns={0 : 'contig', 1 : 'length'})
    data.add_setting('total', df.length.sum())
    data.add_setting('mean', df.iloc[:,1].mean())
    data.add_setting('median', df.iloc[:,1].median())
    data.add_setting('contigs', df.iloc[:,1].count())                                   

 
    
    
    #function to calculate N50
    def N(ratio):
        #break down ratio
        ratio = ratio/100
        #memorizing contig positions
        #start at 0 for single fasta
        mem_pos = [0, (len(df)-1)//2+1]
        #memorizing hits that wer too high & hits that were too low
        too_high = [len(df)-1]
        too_low = [0]
    
    # =============================================================================
    #     analyze the position and the lower position, are we at the N50?
    #     approach: 
    #               if higher jump in middle between current and last lower value
    #               if lower jump in middle between current and last higher value
    # =============================================================================
    
        while 1:

            # N50 at exact position
            if df.iloc[:mem_pos[-1],:].length.sum() >= data.settings['total']*ratio and df.iloc[:mem_pos[-1]-1,:].length.sum() < data.settings['total']*ratio \
            or df.iloc[:mem_pos[-1],:].length.sum() >= data.settings['total']*ratio and len(df) == 1:
                #found the N50!
                if len(df) > 1:
                    N50 = df.iloc[mem_pos[-1],1]
                    break
                #exception for single fasta
                else:
                    N50 = df.iloc[0,1]
                    mem_pos=[0]
                    break
    
            #N50 already passed
            elif df.iloc[:mem_pos[-1],:].length.sum() > data.settings['total']*ratio:
                too_high.append(mem_pos[-1])
                if mem_pos[-1] != mem_pos[-2]:
                    mem_pos.append( mem_pos[-1] - (abs(mem_pos[-1]-too_low[-1])//2))
    
            #N50 not reached yet
            elif df.iloc[:mem_pos[-1],:].length.sum() < data.settings['total']*ratio:
                too_low.append(mem_pos[-1])
                if mem_pos[-1] != mem_pos[-2]:
                    mem_pos.append( mem_pos[-1] + (abs(mem_pos[-1]-too_high[-1])//2))
                #N50 is on last contig
                elif mem_pos[-1]+2 == len(df):
                        N50 = df.iloc[-1,1]
                        break
                #stagnating one lower than N50
                else:
                    mem_pos.append( mem_pos[-1] + 1)
            else:
                print('ERROR: no exception-handling for this type of data!!')
                break

        return N50





    data.add_setting('time3', time.time())

    #create output table
    results = [['contigs:', data.settings['contigs']],
               ['mean:', np.round(data.settings['mean'])],
               ['median:', np.round(data.settings['median'])],
               ['total_bp:', data.settings['total']],
               ['N10:', N(10)],
               ['N50:', N(50)],
               ['N90:', N(90)],
               ['G/C:', '{:.1f}%'.format(GC_count/data.settings['total']*100)]]               

    #small function for printing pandas dataframes
    def pd_print(df):
        df = df.reindex()
        temp = df.iloc[0,0]
        df = df.set_index(df.iloc[:,0])
        df.columns = df.iloc[0,:]
        df = df.reindex(df.index.drop(temp))
        del df[temp]
        del df.index.name
        print(df)

    #print output
    print()
    print('reading_time:', '\t',  '{:.1f}'.format(data.settings['time2']-data.settings['time1']))
    print('total_time:  ', '\t',  '{:.1f}'.format(data.settings['time3']-data.settings['time1']))
    print()
    print()
    print()
    print('results:')
    print()
    pd_print(pd.DataFrame(results))
    print()

    sys.exit()

if __name__ == "__main__":
    main(sys.argv[1:])
