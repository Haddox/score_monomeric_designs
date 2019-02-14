"""
Python script for computing structural metrics relating to the sequence and
ABEGO types
"""

import sys
import Bio.PDB
from numpy import pi, average
penalties = {}
aas = 'ACDEFGHIKLMNPQRSTVWY'
script_dir = path.dirname(__file__)
with open(os.path.join(script_dir, 'penalty_table')) as file:
    for line in file.readlines():
        abego_aa = tuple(line.split()[0:2])
        penalties[abego_aa] = float(line.split()[-2])

def phi_psi_omega_to_abego(phi, psi, omega):
    if psi == None: return 'O'
    if omega == None: omega = 180
    if phi == None: phi=90
    phi = 180 * phi / pi
    psi = 180 * psi / pi
    omega = 180 * omega / pi




    if abs(omega) < 90:
        return 'O'
    elif phi > 0:
        if -100.0 <= psi < 100:
            return 'G'
        else:
            return 'E'
    else:
        if -75.0 <= psi < 50:
            return 'A'
        else:
            return 'B'
    return 'X'

def penalty_to_numeral(value):
    cutoffs=[-3,  -1 , -0.5  , 0  ,  0.2 ,  0.6  , 1.0 ,  1.5  ]
    return '%s' % (len([x for x in cutoffs if x < value]) + 1)

def abego_string(phi_psi_omega):
    out = ''
    for x in phi_psi_omega:
        out += phi_psi_omega_to_abego(x[0], x[1], x[2])
    return out

def main( files ):
    output_pdb = True
    long_output = True
    verbose=False
    short_penalty=False
    if '--verbose' in files:
        verbose=True
    if '--nopdb' in files:
        files.remove('--nopdb')
        output_pdb = False
    if '--short' in files:
        files.remove('--short')
        long_output = False
    if '--shortpenalty' in files:
        files.remove('--shortpenalty')
        short_penalty=True
        long_output=False

    output = ""
    for (filenum, fn) in enumerate(files):
        phi_psi = []
        nres=[]
        seqs = []
        structure = Bio.PDB.PDBParser(QUIET=True).get_structure(fn,fn)
        for chain in structure:
            polypeptides = Bio.PDB.PPBuilder().build_peptides(chain)
            for polypeptide in polypeptides:
                nres.append(len(polypeptide))
                seqs.append(polypeptide.get_sequence())
                if len(nres) > 1:
                    nres[-1] = nres[-1] + nres[-2]
                phi_psi += polypeptide.get_phi_psi_list()
        residues = [res for res in structure.get_residues()]
        phi_psi_omega = []
        #print fn, len(phi_psi), len(residues), nres
        for i in range(len(residues)-1):
            if i +1 in nres:
                omega = None
            else:
                a1 = residues[i]['CA'].get_vector()
                a2 = residues[i]['C'].get_vector()
                a3 = residues[i+1]['N'].get_vector()
                a4 = residues[i+1]['CA'].get_vector()
                omega = Bio.PDB.calc_dihedral(a1,a2,a3,a4)
            phi_psi_omega.append((phi_psi[i][0], phi_psi[i][1], omega))
        phi_psi_omega.append((phi_psi[-1][0], phi_psi[-1][1], None))

        temp_abego_string = abego_string(phi_psi_omega)
        my_seq_string = ''
        my_abego_string = ''
        #print abego_string
        #print
        for s in seqs:
            my_seq_string += str(s)
            my_abego_string += temp_abego_string[0:len(str(s))]
            if str(s) != str(seqs[-1]):
                my_seq_string += ' '
                my_abego_string += ' '
                temp_abego_string = temp_abego_string[len(str(s)):]
                #print
                #print abego_string




        my_penalty_string = '55'
        scores = [0,0]
        best_aa_string = my_seq_string[0:2]
        best_aa_scores = '55'
        for i in range(2,len(my_abego_string)-2):
            if ' ' not in my_abego_string[i-2:i+3]:
                if (my_abego_string[i-1:i+2], my_seq_string[i]) in penalties:
                    term = penalties[(my_abego_string[i-1:i+2], my_seq_string[i])]
                    aa_penalties=[(penalties[my_abego_string[i-1:i+2], aa], aa) for aa in aas]
                    aa_penalties.sort()
                    #print aa_penalties
                    best_aa = aa_penalties[-1][-1]
                    best_aa_string += best_aa
                    best_aa_scores += penalty_to_numeral(  penalties[(my_abego_string[i-1:i+2], best_aa)] )
                    scores.append(-10*term)
                    my_penalty_string += penalty_to_numeral(term)

            elif my_abego_string[i] != ' ':
                scores.append(0)
                my_penalty_string += '5'
                best_aa_string += my_seq_string[i]
                best_aa_scores += '5'
            else:
                my_penalty_string += ' '
                best_aa_string += ' '
                best_aa_scores += ' '


        scores += [0,0]
        if verbose:
            for i in range(len(my_seq_string)):
                output += str(i+1) + ' ' + my_seq_string[i] + ' ' + str(scores[i] / -10.0) + ' '
            output += '\n'

        if long_output:
            output += fn + ' ' + my_seq_string + '\n'
            output += fn + ' ' + my_abego_string + ' '
            output += str( average([x / -10.0 for x in scores]) ) + ' '
            output += str( average([min([x / -10.0, 0]) for x in scores]) ) + '\n'
            output += fn + ' ' + my_penalty_string + '55\n'
            output += 'opt' + (' ' * (len(fn) - 3))  + ' ' + best_aa_string + my_seq_string[-2:] + '\n' #\n\n'
            output += 'opt' + (' ' * (len(fn) - 3))  + ' ' + best_aa_scores + '55 \n\n'
        elif short_penalty:
            output += fn + ' ' + my_seq_string + ' ' + my_abego_string + ' ' + my_penalty_string + '55 '
            output += str( average([x / -10.0 for x in scores]) ) + ' '
            output += str( average([min([x / -10.0, 0]) for x in scores]) ) + '\n'



        else:
            output += fn + ' ' + my_seq_string + ' ' + my_abego_string + ' '
            output += str( average([x / -10.0 for x in scores]) ) + ' '
            output += str( average([min([x / -10.0, 0]) for x in scores]) ) + '\n'


        #j = -1
        #for i in range(len(my_seq_string)):
            #print my_seq_string[i], my_abego_string[i],
            #if my_seq_string[i] != ' ':
            #    j += 1
            #    print scores[j]
            #else:
            #    print

        if output_pdb:
            with open(fn) as file:
                lines = [x for x in file.readlines() if 'ATOM' in x]
            curres = int(lines[0][22:26])
            curres_index = 0

            with open(fn.split('/')[-1] + '.abego_seq.pdb','w') as file:
                for line in lines:
                    myres = int(line[22:26])
                    if myres != curres:
                        curres = myres
                        curres_index += 1
                    outline = line[0:61]
                    outline += '%5s' % ('%.1f' % scores[curres_index])
                    #print '%5s' % '%.1f' % scores[curres_index], len('%5s' % '%.1f' % scores[curres_index])

                    outline += line[66:]
                    file.write(outline)
        if filenum > 0 and filenum % 100 == 0:
            print >>sys.stderr, 'Processed', filenum, '/', len(files), 'structures.'
    return output.rstrip('\n')

if __name__ == "__main__":
    print main( sys.argv[1:] )
