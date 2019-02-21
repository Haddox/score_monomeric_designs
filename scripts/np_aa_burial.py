#!/usr/bin/python
import sys
import commands
import time
from multiprocessing import Pool
import math

def burial_info(x):
    fn, np_total, aa_dict = x

    # TODO: Right now the below line uses RosettaScripts to run an XML. It has
    # external dependencies. We should change this so that it works internally,
    # and use the pyrosetta.distributed module if possible.
    out = commands.getoutput('/software/rosetta/versions/v2019.01-dev60566/bin/rosetta_scripts.hdf5.linuxgccrelease -parser:protocol ./scripts/residue_burial_nodes.xml -s %s -nooutput' % fn).split('\n')

    out = [x.strip() for x in out]

    start_nonpolar = out.index('protocols.rosetta_scripts.ParsedProtocol: =======================BEGIN FILTER nonpolar=======================')
    end_nonpolar = out.index('protocols.rosetta_scripts.ParsedProtocol: =======================END FILTER nonpolar=======================')

    sasa = [float(x.split()[-1]) for x in filter(lambda y: 'HYDROPHOBIC SASA' in y, out[start_nonpolar:end_nonpolar])]
    res = [x.split()[1][0:3] for x in filter(lambda y: 'HYDROPHOBIC SASA' in y, out[start_nonpolar:end_nonpolar])]

    buried_np_AFILMVWY=0
    buried_np=0

    exposed_np_AFILMVWY=0
    exposed_np=0

    seq = ''

    for i in range(len(res)):
        buried_np += np_total[res[i]] - sasa[i]
        exposed_np += sasa[i]
        if res[i] in ['PHE','ILE','LEU','MET','VAL','TRP','TYR','ALA']:
            buried_np_AFILMVWY += np_total[res[i]] - sasa[i]
            exposed_np_AFILMVWY += sasa[i]
        seq += aa_dict[res[i]]

    return fn, seq, buried_np, exposed_np, buried_np_AFILMVWY, exposed_np_AFILMVWY

def main( args ):

    data = """
    ALA 111
    CYS 125
    ASP 92
    GLU 124
    PHE 210
    GLY 92
    HIS 155
    ILE 182
    LYS 182
    LEU 183
    MET 208
    ASN 90
    PRO 171
    GLN 118
    ARG 141
    SER 96
    THR 121
    VAL 157
    TRP 222
    TYR 177""".split('\n')[1:]
    np_total = {}

    for line in data:
        aa = line.split()[0]
        np_total[aa] = float(line.split()[1])

    aas = """A ALA
        C CYS
        D ASP
        E GLU
        F PHE
        G GLY
        H HIS
        I ILE
        K LYS
        L LEU
        M MET
        N ASN
        P PRO
        Q GLN
        R ARG
        S SER
        T THR
        V VAL
        W TRP
        Y TYR""".split('\n')
    aa_dict = {}
    for aa in aas:
        aa_dict[aa.split()[1]] = aa.split()[0]


    files = args
    output = ''
    if '--noheader' in files:
        files.remove('--noheader')
        output = ''
    else:
        output = 'description   sequence   buried_np   exposed_np   buried_np_AFILMVWY   exposed_np_AFILMVWY\n'

    for fn in files:
        (fn_out, seq, buried_np, exposed_np, buried_np_AFILMVWY, exposed_np_AFILMVWY) = burial_info([fn, np_total, aa_dict])
        output += '%s   %s   %s   %s   %s   %s\n' % (fn_out.split('/')[-1], seq, buried_np, exposed_np, buried_np_AFILMVWY, exposed_np_AFILMVWY)

    return output.rstrip( '\n' )

if __name__ == "__main__":
    print main( sys.argv[1:] )
