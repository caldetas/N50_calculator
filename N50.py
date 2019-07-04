#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 29 14:34:57 2018

@author: hannes
"""

import sys
import getopt
import time


#instructions
form='N50.py\n\t--minc *minimum contig length, opt*\n\t <flat>'

time1 = time.time()

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
            temp=''
        else:
            temp += line

    #last contig
    l1[-1].append(len(temp))

    temp=''
    t=0

    ass.close()
   
    time2 = time.time()
    print()
    print('reading_time:', np.round(time2-time1))
    print()
    print()

    df = pd.DataFrame(l1)

    df = df.iloc[:,:2]
    df = df.sort_values(by=[1])
    df = df.reset_index(drop=True)
    
    df = df.loc[df[1] >= minc]
    

    df = df.rename(columns={0 : 'contig', 1 : 'length'})#, 2 : 'subtotal'})
    total = df.length.sum()

    
    #memorizing contig positions, jump always in middle of ordered df
    #start at 0
    mem_pos = [0, (len(df)-1)//2+1]

# =============================================================================
#     analyze last two positions, are they in N50 region?
#     approach: jump to middle value
#     if not at N50: add half of position-difference
#     if over N50: subtract half of position-difference
#     check if reached: is lower neighbour already in list?    
# =============================================================================

    while 1:
#        print(mem_pos)
#        print(df.iloc[:1,:])
#        print(mem_pos[-10:])
#        print(df.iloc[mem_pos[-1],:])

        if df.iloc[:mem_pos[-1]+1,:].length.sum() >= total*0.5 and mem_pos[-1]-1 in mem_pos\
        or df.iloc[:mem_pos[-1]+1,:].length.sum() >= total*0.5 and mem_pos[-1]+1 in mem_pos:
#            print('found')
            #found the N50!
            if len(mem_pos) > 2:
                N50 = df.iloc[mem_pos[-1],1]
                break
            #exception for single fasta
            else:
                if df.iloc[:mem_pos[-2]+1,:].length.sum() >= total*0.5:
                    N50 = df.iloc[0,1]
                    mem_pos=mem_pos[:1]
                    break

        #N50 already reached
        elif df.iloc[:mem_pos[-1]+1,:].length.sum() > total*0.5:
#            print('higher')
            #print(df.iloc[:mem_pos[-1]+1,:].length.sum(), total*.5)
            mem_pos.append( mem_pos[-1] - (abs(mem_pos[-1]-mem_pos[-2])//2))
            

        #N50 not reached yet
        elif df.iloc[:mem_pos[-1]+1,:].length.sum() < total*0.5:
#            print('lower')
            if mem_pos[-1] != mem_pos[-2]:
                mem_pos.append( mem_pos[-1] + (abs(mem_pos[-1]-mem_pos[-2])//2))
            #stagnating on lower values
            elif df.iloc[:mem_pos[-1]+2,:].length.sum() >= total*0.5:
#                print('exeption2')
                N50 = df.iloc[mem_pos[-1]+2,1]
                break
            #stagnating on lower values
            elif df.iloc[:mem_pos[-1]+3,:].length.sum() >= total*0.5:
#                print('exeption3')
                N50 = df.iloc[mem_pos[-1]+3,1]
                break

        #exception: exact N50 nr gets hit from a not neighbouring position
        elif df.iloc[:mem_pos[-1]+1,:].length.sum() == total*0.5 and df.iloc[:mem_pos[-1],:].length.sum() < total*0.5:
#            print('exeption1')
            #found the N50!)
            N50 = df.iloc[mem_pos[-1],1]
            #exit the while loop
            break
        else:
            print('ERROR: no exception-handling for this type of data!!')
            break

    mean = df.iloc[:,1].mean()
    median = df.iloc[:,1].median()
    contigs =  df.iloc[:,1].count()
    time3 = time.time()
    
    results = [['total_time:', np.round(time3-time1)],
                ['contigs:', contigs],
               ['mean:', np.round(mean)],
               ['median:', np.round(median)],
               ['total_bp:', total],
               ['N50:', N50]]
    results_pd = pd.DataFrame(results, dtype=int)
    
    results_pd = results_pd.set_index([0])
    results_pd.columns = results_pd.iloc[0]
    results_pd = results_pd.reindex(results_pd.index.drop('total_time:'))
    del results_pd.index.name



    for i in range(len(mem_pos)):
        mem_pos[i] += 1
  
    print('N50-search history on contigs:')
    print()
    print(mem_pos)
    print()
    print()
    print('results:')
    print()
    print(results_pd)
    print()

    



        
        
    sys.exit()

if __name__ == "__main__":
    main(sys.argv[1:])
