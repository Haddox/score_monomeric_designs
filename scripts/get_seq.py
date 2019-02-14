#!/usr/bin/python
import commands
import sys
import os

def main( args ):
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

    seqs = {}

    chainA = False
    chainB = False
    model1 = False
    skipcompare = False
    files = args
    if 'A' in files and not os.path.isfile('A'):
        chainA = True
        files.remove('A')

    if 'B' in files and not os.path.isfile('B'):
        chainB = True
        files.remove('B')

    if '1' in files and not os.path.isfile('1'):
        model1 = True
        files.remove('1')

    if '--seqonly' in files:
        skipcompare = True
        files.remove('--seqonly')

    for (i, filename) in enumerate(sorted(files)):
        with open(filename) as file:
            lines = file.readlines()
            lastline = 0
            for lastline in range(len(lines)):
                if len(lines[lastline]) > 6:
                    if 'ENDMDL' in lines[lastline][0:6]: break

            lines = lines[0:lastline]

            lines = filter(lambda x: ' CA ' in x and 'ATOM' in x, lines)
            if chainA: lines = filter(lambda x: ' A ' in x, lines)
            if chainB: lines = filter(lambda x: ' B ' in x, lines)
        seq = ''
        for line in lines:
            seq += aa_dict[line.split()[3]]
        seqs[filename] = seq
        if i > 0 and i % 1000 == 0: print >>sys.stderr, 'Read', i, 'of', len(files)

    n = max([len(x) for x in seqs])

    output = ""
    for seq in seqs:
        output += seq + (' ' * (n-len(seq))) + ' ' + seqs[seq] + '\n'
#print (' '* n),

    if not skipcompare:
        diffs = ''
        for i in range(len(seqs[seq])):
            diffs += ' '
            for seqb in seqs:
                #try:
                    if seqs[seqb][i] != seqs[seq][i]:
                        diffs = diffs[0:-1] + '*'
                #except:
                #    print seqb, seqs[seqb], 'wtf'
        output += diffs
    return output.rstrip( '\n' )

if __name__ == "__main__":
    print main( sys.argv[1:] )
