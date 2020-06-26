#!/usr/bin/python
from numpy import *
from matplotlib import *
import matplotlib
matplotlib.use('Agg')
from pylab import *
import sys
import os

pwd = os.path.abspath('.')
files = sys.argv[1:]
plotme = True
if '--noplot' in files:
    files.remove('--noplot') 
    plotme = False

if files == []: files = ['start_frags.fsc.100.9mers']

for fn in files:
    if 'frag_qual.dat' not in fn: continue
    average_all_frags = []

    #print os.path.abspath('.')
    directory = fn.replace('frag_qual.dat','')
    if directory != '': os.chdir(directory)

    xs, ys=[],[]
    best = {}
    if os.path.isfile('frag_qual.summary'):
        with open('frag_qual.summary') as file:
            lines = file.readlines()
        print lines[0].strip()
        os.chdir(pwd)
        continue
    with open('frag_qual.dat') as file:
        lines = file.readlines()
    for line in lines:
        xs.append(int(line.split()[1]))
        ys.append(float(line.split()[3]))
        if xs[-1] not in best:
            best[xs[-1]] = []
        best[xs[-1]].append(ys[-1])
        average_all_frags.append(ys[-1])
    for i in best:
        best[i].sort()
    lowest = [best[x][0] for x in best]
    lowest.sort()
    lowest.reverse()
    if plotme and not os.path.isfile('frag_qual.png'):
        Figure()
        plot(xs,ys,'.')
        ylim(0,5)
        savefig('frag_qual.png')
        close()

    with open('00001.fasta') as file:
        seq = file.readlines()[-1].strip()


    print fn, seq, lowest[0], lowest[1], lowest[2], lowest[3], lowest[4], lowest[5], 'sum6', sum(lowest[0:6]), 'sum', sum(lowest), average(average_all_frags)
    with open('frag_qual.summary','w') as file:
        file.write('%s ' % fn)
        file.write('%s ' % seq)
        file.write('%s ' % lowest[0])
        file.write('%s ' % lowest[1])
        file.write('%s ' % lowest[2])
        file.write('%s ' % lowest[3])
        file.write('%s ' % lowest[4])
        file.write('%s ' % lowest[5])
        file.write('sum6 ')
        file.write('%s ' % sum(lowest[0:6]))
        file.write('sum ')
        file.write('%s ' % sum(lowest))
        file.write('%s\n' % average(average_all_frags))

    os.chdir(pwd)
