#!/usr/bin/python
import sys
import commands
from numpy import *
from numpy.linalg import norm
import get_seq

def main( args ):
    outstr = ''
    frag_picker_dir = '/software/frapicker/'
    for filenum, fn in enumerate(args):
        seq = get_seq.main([fn]).split('\n')[0].split()[-1]
        COMs = []
        with open(fn) as file:
            lines = file.readlines()
        lines = filter(lambda x: ' CA ' in x, lines)
        dssp_string = commands.getoutput('{0}/Rosetta/tools/fragment_tools/pdb2vall/structure_profile_scripts/dssp2threestateSS.pl {1}'.format(frag_picker_dir, fn)).split('\n')[-1].strip()
        if 'E' in dssp_string:
            outstr += '%s %s %s\n' % (fn, seq, 0)
            continue
        #print dssp_string
        helices = [[]]
        cur_helix = 0
        helix_indices = [x for x in range(len(dssp_string)) if dssp_string[x] == 'H']
        helices[0].append(helix_indices[0])
        for i in range(1,len(helix_indices)):
            if helix_indices[i] != helix_indices[i-1] + 1:
                helices.append([])
            #if i == 31 and len(helices) == 2: helices.append([])
            helices[-1].append(helix_indices[i])
        if len(helices) != 3:
            outstr += '%s %s %s\n' % (fn, seq, 0)
            continue
        for helix in helices:
            sumcoords = array([0.0, 0.0, 0.0])
            for res in helix:
                sumcoords += array([float(x) for x in lines[res].split()[6:9]])

            COMs.append( sumcoords / len(helix))
        #print COMs
        v1 = COMs[1] - COMs[0]
        v2 = COMs[2] - COMs[0]
        cross1 = cross(v1,v2)
        #print cross1, norm(cross1)
        cross1 = cross1 / norm(cross1)


        helix1_0 = array([float(x) for x in lines[helices[0][0]].split()[6:9]])
        helix1_1 = array([float(x) for x in lines[helices[0][-1]].split()[6:9]])
        v3 = helix1_1 - helix1_0
        v3 = v3 / norm(v3)
        #print cross1, v3
        #seq = commands.getoutput('get_seq.py %s' % fn).split('\n')[-2].strip()
        outstr += '%s %s %s %s %s %s %s\n' % (fn, seq, dot(cross1, v3), len(helices[0]), len(helices[1]), len(helices[2]), dssp_string)
        if filenum > 0 and filenum % 100 == 0:
            print >>sys.stderr, 'Processed', filenum, '/', len(args), 'structures.'
    return outstr.strip()

if __name__ == "__main__":
    print main( sys.argv[1:] ),
