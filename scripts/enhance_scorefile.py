"""
Python script for computing a variety of structural metrics
"""

import numpy as np
import pandas as pd
import sys
import os
import json

import np_aa_burial
# import commands
import subprocess
import get_seq
import abego_seq_agreement
import core_clusters
import helix_handedness

score_dict = pd.read_csv(sys.argv[1]) #, header=1, sep = '\s+')
pdbs = ['%s.pdb' % x for x in score_dict['description'].values]
pdbs = ' '.join(pdbs)

if os.path.isfile('%s.sequences' % sys.argv[1]):
    print >> sys.stderr, 'Loading sequences from %s.sequences' % sys.argv[1], sys.stderr
    with open('%s.sequences' % sys.argv[1]) as file:
        lines = file.readlines()
else:
    print >> sys.stderr, 'Collecting sequences...', sys.stderr
    lines = get_seq.main( [ '--seqonly' ] + pdbs.split() ).split( '\n' )
    with open('%s.sequences' % sys.argv[1],'w') as file:
        for line in lines:
            file.write('%s\n' % line)
name_to_seq = {}
for line in lines:
    name = line.split()[0][0:-4]
    seq = line.split()[-1].strip()
    name_to_seq[name] = seq
print >> sys.stderr, (name_to_seq), sys.stderr
print >> sys.stderr, (score_dict['description']), sys.stderr
score_dict['sequence'] = score_dict['description'].map(lambda x: name_to_seq[x])
score_dict['n_res'] = score_dict['sequence'].map(len)


#dssp
if os.path.isfile('%s.dssp' % sys.argv[1]):
    print >> sys.stderr, 'Loading DSSP from %s.dssp' % sys.argv[1], sys.stderr
    with open('%s.dssp' % sys.argv[1]) as file:
        lines = file.readlines()
    name_to_dssp = {}
    for line in lines:
        name = line.split()[0]
        dssp = line.split()[-1].strip()
        name_to_dssp[name] = dssp
else:
    print >> sys.stderr, 'Collecting DSSP...', sys.stderr
    name_to_dssp = {}
    for (i,pdb) in enumerate(pdbs.split()):
        # TODO fixup dssp calc to not write extra .dssp files
        # Copy pdb file to working directory
        # Run dssp on the tempfile
        cmd = ['/work/robetta/workspace/labFragPicker_DO_NOT_REMOVE/Rosetta/tools/fragment_tools/pdb2vall/structure_profile_scripts/dssp2threestateSS.pl', pdb]
        process = subprocess.Popen(cmd, stdout=subprocess.PIPE)
        out, err = process.communicate()
        name_to_dssp[pdb[0:-4]] = ''.join( out.split('\n')[1:] )
        #name_to_dssp[pdb[0:-4]] = ''.join( commands.getoutput('/work/robetta/workspace/labFragPicker_DO_NOT_REMOVE/Rosetta/tools/fragment_tools/pdb2vall/structure_profile_scripts/dssp2threestateSS.pl %s' % pdb).split('\n')[1:] )
        #print >> sys.stderr, pdb, name_to_dssp[pdb[0:-4]]
        if i % 100 == 0:
            print >> sys.stderr, 'Obtained DSSP for structure %d of %d' % ( i+1, len(pdbs.split()) ), sys.stderr
    subprocess.call(['rm', '*_0001.pdb.dssp'])
    with open('%s.dssp' % sys.argv[1],'w') as file:
        for pdb in pdbs.split():
            file.write('%s %s\n' % (pdb[0:-4], name_to_dssp[pdb[0:-4]]))
score_dict['dssp'] = score_dict['description'].map(lambda x: name_to_dssp[x])
score_dict['nres_helix'] = score_dict['dssp'].map(lambda x: len([r for r in x if r == 'H']))
score_dict['nres_sheet'] = score_dict['dssp'].map(lambda x: len([r for r in x if r == 'E']))
score_dict['nres_loop'] = score_dict['dssp'].map(lambda x: len([r for r in x if r == 'L']))
score_dict['frac_helix'] = map(lambda x, y: float(x) / float(y), score_dict['nres_helix'], score_dict['n_res'])
score_dict['frac_sheet'] = map(lambda x, y: float(x) / float(y), score_dict['nres_sheet'], score_dict['n_res'])
score_dict['frac_loop'] = map(lambda x, y: float(x) / float(y), score_dict['nres_loop'], score_dict['n_res'])
score_dict['hbond_sr_bb_per_helix'] = map(lambda x, y:  y/max([1,y]) *      (float(x) / max([1,float(y)])), score_dict['hbond_sr_bb'], score_dict['nres_helix'])
score_dict['hbond_lr_bb_per_sheet'] = map(lambda x, y:  y/max([1,y]) *      (float(x) / max([1,float(y)])), score_dict['hbond_lr_bb'], score_dict['nres_sheet'])



def helix_terms(seq, dssp):
    helices=[[-10]]
    for i in range(len(seq)):
        if dssp[i] == 'H':
            if i == helices[-1][-1] + 1:
                helices[-1].append(i)
            else:
                helices.append([i])
    T1=[]
    Tminus1=[]
    for H in helices:
        if len(H) >= 5:
            T1 += H[0:3]
            Tminus1 += H[-3:]
    T1aa = [seq[x] for x in T1]
    Tminus1aa = [seq[x] for x in Tminus1]
    def netQ(seq):
        qs=dict(K=1,
                R=1,
                D=-1,
                E=-1,
                H=0.5,)
        return sum([qs[x] for x in seq if x in qs])
    def absQ(seq):
        qs=dict(K=1,
                R=1,
                D=1,
                E=1,
                H=0.5,)
        return sum([qs[x] for x in seq if x in qs])


    return [netQ(T1aa), netQ(Tminus1aa), netQ(Tminus1aa) - netQ(T1aa), absQ(T1aa), absQ(Tminus1aa), absQ(T1aa) + absQ(Tminus1aa)]
helix_netq=np.array([helix_terms(seq,dssp) for seq,dssp in zip(score_dict['sequence'], score_dict['dssp'])])
print >> sys.stderr, np.shape(helix_netq), sys.stderr
score_dict['T1_netq'] = helix_netq[:,0]
score_dict['Tminus1_netq'] = helix_netq[:,1]
score_dict['Tend_netq'] = helix_netq[:,2]
score_dict['T1_absq'] = helix_netq[:,3]
score_dict['Tminus1_absq'] = helix_netq[:,4]
score_dict['Tend_absq'] = helix_netq[:,5]


#gc_content
if os.path.isfile('%s.gc_content' % sys.argv[1]):
    print >> sys.stderr, 'Loading GC content data from %s.gc_content...' % sys.argv[1], sys.stderr
    seq_to_gc_content = {}
    with open('%s.gc_content' % sys.argv[1]) as file:
        gclines = file.readlines()
        for line in gclines:
            seq = line.split()[0]
            seq_to_gc_content[seq] = float(line.split()[1])

    score_dict['GC_content'] = score_dict['sequence'].map(lambda x: seq_to_gc_content[x])
else:
    print >> sys.stderr, '%s.gc_content not found, skipping GC content' % sys.argv[1], sys.stderr


#fragments
if os.path.isfile('%s.frag_quals' % sys.argv[1]):
    print >> sys.stderr, 'Loading fragment quality data from %s.frag_quals...' % sys.argv[1], sys.stderr
    seq_to_frag_line = {}
    with open('%s.frag_quals' % sys.argv[1]) as file:
        fraglines = file.readlines()
        for line in fraglines:
            seq = line.split()[1]
            if seq in seq_to_frag_line:
                #TODO refactor fragment quality assesment system to use conformation dependent
                #identifier. In reality this script should be calling make-fragments on a per-pdb
                #basis to generate fragment scores instead of relying on an external
                #make fragments call.
                raise ValueError("Redundant sequence found in frag_quals file.")
            seq_to_frag_line[seq] = line

    #print >> sys.stderr, len(score_dict)
    score_dict['found_sequence_for_frags'] = score_dict['sequence'].map(lambda x: x in seq_to_frag_line)
    score_dict = score_dict[score_dict.found_sequence_for_frags == True]
    score_dict = score_dict.drop('found_sequence_for_frags', 1)
    #print >> sys.stderr, len(score_dict)

    score_dict['worstfrag'] = score_dict['sequence'].map(lambda x: float(seq_to_frag_line[x].split()[2]))
    score_dict['worst6frags'] = score_dict['sequence'].map(lambda x: float(seq_to_frag_line[x].split()[-4]))
    score_dict['sum_best_frags'] = score_dict['sequence'].map(lambda x: float(seq_to_frag_line[x].split()[-2]))
    score_dict['avg_all_frags'] = score_dict['sequence'].map(lambda x: float(seq_to_frag_line[x].split()[-1]))
    n_frags = float(score_dict['n_res']) - 8
    score_dict['avg_best_frag'] = float(score_dict['sum_best_frags']) / n_frags
else:
    print >> sys.stderr, '%s.frag_quals not found, skipping fragment quality' % sys.argv[1], sys.stderr

#abego terms
seq_to_abego_terms = {}
if os.path.isfile('%s.abego_terms' % sys.argv[1]):
    print >> sys.stderr, 'Loading abego terms from %s.abego_terms...' % sys.argv[1], sys.stderr
    with open('%s.abego_terms' % sys.argv[1]) as file:
        lines = file.readlines()
else:
    print >> sys.stderr, 'Calculating abego term penalties...', sys.stderr
    #print >> sys.stderr, '/work/grocklin/gabe_scripts/abego_seq_agreement.py --nopdb --short %s > %s.abego_terms' % (pdbs,sys.argv[1])
    #lines = commands.getoutput('/work/grocklin/gabe_scripts/abego_seq_agreement.py --nopdb --short %s > %s.abego_terms' % (pdbs,sys.argv[1])).split('\n')
    lines = abego_seq_agreement.main( [ '--nopdb', '--short' ] + pdbs.split() ).split( '\n' )
    #print >> sys.stderr, lines
    with open('%s.abego_terms' % sys.argv[1],'w') as file:
        for line in lines:
            file.write('%s\n' % line)
for line in lines:
    seq = line.split()[1]
    seq_to_abego_terms[seq] = [float(x) for x in line.split()[-2:]]

score_dict['abego_res_profile'] = score_dict['sequence'].map(
    lambda x: float(seq_to_abego_terms[x][-2]) if x in seq_to_abego_terms.keys() else "KeyError: {0}".format(x)
)
score_dict['abego_res_profile_penalty'] = score_dict['sequence'].map(
    lambda x: float(seq_to_abego_terms[x][-1]) if x in seq_to_abego_terms.keys() else "KeyError: {0}".format(x)
)

#helix handedness
# seq_to_3helix_handedness = {}
# if os.path.isfile('%s.3helix_handedness' % sys.argv[1]):
#     print >> sys.stderr, 'Loading 3-helix bundle handedness from %s.3helix_handedness...' % sys.argv[1]
#     with open('%s.3helix_handedness' % sys.argv[1]) as file:
#         lines = file.readlines()
# else:
#     print >> sys.stderr, 'Calculating 3-helix bundle handedness...'
#     #print >> sys.stderr, '/work/grocklin/gabe_scripts/abego_seq_agreement.py --nopdb --short %s > %s.abego_terms' % (pdbs,sys.argv[1])
#     #lines = commands.getoutput('/work/grocklin/gabe_scripts/abego_seq_agreement.py --nopdb --short %s > %s.abego_terms' % (pdbs,sys.argv[1])).split('\n')
#     lines = helix_handedness.main( pdbs.split() ).split( '\n' )
#     #print >> sys.stderr, lines
#     with open('%s.3helix_handedness' % sys.argv[1],'w') as file:
#         for line in lines:
#             file.write('%s\n' % line)
# for line in lines:
#     seq = line.split()[1]
#     seq_to_3helix_handedness[seq] = float(line.split()[2])
# score_dict['3helix_handedness'] = score_dict['sequence'].map(lambda x: float(seq_to_3helix_handedness[x]))


#core cluster terms
seq_to_core_cluster_terms = {}
if os.path.isfile('%s.core_clusters' % sys.argv[1]):
    print >> sys.stderr, 'Loading core cluster terms from %s.core_clusters...' % sys.argv[1], sys.stderr
    with open('%s.core_clusters' % sys.argv[1]) as file:
        lines = file.readlines()
else:
    print >> sys.stderr, 'Calculating core cluster terms...', sys.stderr
    #print >> sys.stderr, pdbs
    #lines = commands.getoutput('/work/grocklin/gabe_scripts/core_clusters.py %s > %s.core_clusters' % (pdbs, sys.argv[1])).split('\n')
    lines = core_clusters.main( pdbs.split() ).split( '\n' )
    with open('%s.core_clusters' % sys.argv[1],'w') as file:
        for line in lines:
            file.write('%s\n' % line)
for line in lines:
    seq = line.split()[1]
    seq_to_core_cluster_terms[seq] = line.split()[2:6]

score_dict['largest_hphob_cluster'] = score_dict['sequence'].map(lambda x: int(seq_to_core_cluster_terms[x][0]))
score_dict['n_hphob_clusters'] = score_dict['sequence'].map(lambda x: int(seq_to_core_cluster_terms[x][1]))
score_dict['hphob_sc_contacts'] = score_dict['sequence'].map(lambda x: int(seq_to_core_cluster_terms[x][2]))
score_dict['hphob_sc_degree'] = score_dict['sequence'].map(lambda x: float(seq_to_core_cluster_terms[x][3]))

#forward folding terms
seq_to_ffmetric = {}
if os.path.isfile('%s.ffmetric_detailed' % sys.argv[1]):
    print >> sys.stderr, 'Loading forward folding metric from %s.ffmetric_detailed...' % sys.argv[1], sys.stderr
    ffmetric_detailed = pd.read_csv('%s.ffmetric_detailed' % sys.argv[1],delim_whitespace=True)
    for goodcol in 'ffmetric  ffmetric_withrelax  low_rms_abinit  min_energy_abinit  min_energy_rms_abinit  min_energy_all  min_energy_rms_all egap1A  egap2A  egap3A  egap_relax1A  egap_relax2A  egap_relax3A'.split():
        datadict = ffmetric_detailed.set_index('sequence').to_dict()[goodcol]
        for seq in score_dict['sequence'].values:
            if seq not in datadict:
                datadict[seq] = np.median(ffmetric_detailed[goodcol].values)
        score_dict[goodcol] = score_dict['sequence'].map(lambda x: datadict[x])

elif os.path.isfile('%s.ffmetric' % sys.argv[1]):
    print >> sys.stderr, 'Loading forward folding metric from %s.ffmetric...' % sys.argv[1], sys.stderr
    with open('%s.ffmetric' % sys.argv[1]) as file:
        lines = file.readlines()
    for line in lines:
        seq = line.split()[0]
        seq_to_ffmetric[seq] = float(line.split()[1])
    score_dict['ffmetric'] = score_dict['sequence'].map(lambda x: seq_to_ffmetric[x])
else:
    print >> sys.stderr, '%s.ffmetric not found, skipping forward folding' % sys.argv[1], sys.stderr

#gevorg terms
if os.path.isfile('%s.term_analysis' % sys.argv[1]):
    print >> sys.stderr, 'Loading TERM analysis from %s.term_analysis...' % sys.argv[1], sys.stderr
    description_to_gevorg_metric = {}
    with open('%s.term_analysis' % sys.argv[1]) as file:
        lines = file.readlines()
        for line in lines:
            description = line.split()[0] + '_0001'
            description_to_gevorg_metric[description] = [float(x) for x in line.split()[1:]]
        for description in score_dict.description.values:
            if description not in description_to_gevorg_metric: description_to_gevorg_metric[description] = [0,0,0,0,0,0]
    score_dict['abd50_mean'] = score_dict['description'].map(lambda x: description_to_gevorg_metric[x][0])
    score_dict['abd50_min'] = score_dict['description'].map(lambda x: description_to_gevorg_metric[x][1])
    score_dict['dsc50_mean'] = score_dict['description'].map(lambda x: description_to_gevorg_metric[x][2])
    score_dict['dsc50_min'] = score_dict['description'].map(lambda x: description_to_gevorg_metric[x][3])
    score_dict['ssc50_mean'] = score_dict['description'].map(lambda x: description_to_gevorg_metric[x][4])
    score_dict['ssc50_min'] = score_dict['description'].map(lambda x: description_to_gevorg_metric[x][5])
else:
    print >> sys.stderr, '%s.term_analysis not found, skipping TERM analysis' % sys.argv[1], sys.stderr

#nonlocal analysis
if os.path.isfile('%s.nonlocal' % sys.argv[1]):
    print >> sys.stderr, 'Loading nonlocal analysis from %s.nonlocal...' % sys.argv[1], sys.stderr
    seq_to_nonlocal = {}
    with open('%s.nonlocal' % sys.argv[1]) as file:
        lines = file.readlines()
        for line in lines:
            seq = line.split()[1]
            seq_to_nonlocal[seq] = [float(line.split()[x]) for x in [3,4,6,7]]
        for seq in score_dict.sequence.values:
            if seq not in seq_to_nonlocal: seq_to_nonlocal[seq] = [0,0,0,0]
    score_dict['nonlocal_total_per_res'] = score_dict['sequence'].map(lambda x: seq_to_nonlocal[x][0])
    score_dict['nonlocal_total_fxn'] = score_dict['sequence'].map(lambda x: seq_to_nonlocal[x][1])
    score_dict['nonlocal_fa_atr_per_res'] = score_dict['sequence'].map(lambda x: seq_to_nonlocal[x][2])
    score_dict['nonlocal_fa_atr_fxn'] = score_dict['sequence'].map(lambda x: seq_to_nonlocal[x][3])
else:
    print >> sys.stderr, '%s.nonlocal not found, skipping nonlocal analysis' % sys.argv[1], sys.stderr


#residue-specific burial analysis
if not os.path.isfile('%s.burial_analysis' % sys.argv[1]):
    print >> sys.stderr, 'Calculating burial analysis...', sys.stderr
    lines = np_aa_burial.main( pdbs.split() )
    with open('%s.burial_analysis' % sys.argv[1],'w') as file:
        file.write(lines)

else:
    print >> sys.stderr, 'Loading burial analysis from %s.burial_analysis...' % sys.argv[1], sys.stderr
burial_df = pd.read_csv('%s.burial_analysis' % sys.argv[1],delim_whitespace=True)
score_dict = pd.merge(left=score_dict, right=burial_df[['sequence','buried_np_AFILMVWY','exposed_np_AFILMVWY']], on='sequence',how='left')
score_dict['buried_np_AFILMVWY_per_res'] = float(score_dict.buried_np_AFILMVWY) / float(score_dict.n_res)


#residue scores vs means analysis
# seq_to_score_vs_mean = {}
# if os.path.isfile('%s.score_vs_mean' % sys.argv[1]):
#     print >> sys.stderr, 'Loading residue scores vs means analysis from %s.score_vs_mean...' % sys.argv[1]
# else:
#     print >> sys.stderr, 'Calculating residue scores vs means...'
#     commands.getoutput('/work/grocklin/gabe_scripts/score_vs_mean.py %s --quiet > %s.score_vs_mean' % (pdbs, sys.argv[1])).split('\n')
# with open('%s.score_vs_mean' % sys.argv[1]) as file:
#     lines = file.readlines()
#     for line in lines:
#         seq = line.split()[1]
#         seq_to_score_vs_mean[seq] = [float(line.split()[x]) for x in [2,3,4]]
#     for seq in score_dict.sequence.values:
#         if seq not in seq_to_score_vs_mean: seq_to_score_vs_mean[seq] = [0,0,0]
# score_dict['res_vs_mean'] = score_dict['sequence'].map(lambda x: seq_to_score_vs_mean[x][0])
# score_dict['res_vs_mean_worst'] = score_dict['sequence'].map(lambda x: seq_to_score_vs_mean[x][1])
# score_dict['res_vs_mean_4th'] = score_dict['sequence'].map(lambda x: seq_to_score_vs_mean[x][2])

#net charge
charges = {'H': 0.5, 'R': 1.0, 'K': 1.0, 'E': -1.0, 'D': -1.0}
def seq_to_net_charge(seq):
    out = 0.0
    for s in seq:
        if s in charges:
            out += charges[s]
    return out
score_dict['netcharge'] = score_dict['sequence'].map(seq_to_net_charge)

#n_charged
def seq_to_n_charged(seq):
    return float(len([x for x in seq if x in 'RKDE'])) + (0.5 * float(len([x for x in seq if x == 'H'])))
score_dict['n_charged'] = score_dict['sequence'].map(seq_to_n_charged)


#topology_codes = {'HHH':1, 'HEEH': 2, 'EHEE': 3, 'EEHEE': 4}
#def name_to_topology_code(name):
#    if name.split('_')[0] in topology_codes:
#        return topology_codes[name.split('_')[0]]
#    else:
#        return 'na'
#score_dict['topology_code'] = score_dict['description'].map(name_to_topology_code)
#
#def name_to_rank(name):
#    if len(name.split('_')) < 2: return 'na'
#    if name.split('_')[-2].isdigit():
#        return int(name.split('_')[-2])
#    else:
#        return 'na'
#score_dict['rank'] = score_dict['description'].map(name_to_rank)

#n_hydrophobic
def seq_to_n_hydrophobic(seq):
    return len([x for x in seq if x in 'FIWLVMYA'])
def seq_to_n_hydrophobic_noA(seq):
    return len([x for x in seq if x in 'FIWLVMY'])
score_dict['n_hydrophobic'] = score_dict['sequence'].map(seq_to_n_hydrophobic)
score_dict['n_hydrophobic_noA'] = score_dict['sequence'].map(seq_to_n_hydrophobic_noA)

hydrophobicity_data="""F 100
I 99
W 97
L 97
V 76
M 74
Y 63
C 49
A 41
T 13
H 8
G 0
S -5
Q -10
R -14
K -23
N -28
E -31
P -46
D -55"""
hydrophobicity_score={}
for line in hydrophobicity_data.split('\n'):
    aa = line.split()[0]
    score = float(line.split()[1])
    hydrophobicity_score[aa]=score

def seq_to_hydrophobicity(seq):
    out = 0
    for aa in seq:
        out += hydrophobicity_score[aa]
    return out
score_dict['hydrophobicity'] = score_dict['sequence'].map(seq_to_hydrophobicity)


def seq_to_contig_not_hp(seq):
    nonhp = ''
    for aa in seq:
        if aa in 'FIWLVMY':
            nonhp += ' '
        else:
            nonhp += aa
    return max([len(x) for x in nonhp.split()])
score_dict['contig_not_hp_max'] = score_dict['sequence'].map(seq_to_contig_not_hp)


def seq_to_contig_not_hp_internal(seq):
    nonhp = ''
    for aa in seq:
        if aa in 'FIWLVMY':
            nonhp += ' '
        else:
            nonhp += aa
    contigs = [len(x) for x in nonhp.split()]
    if len(contigs) >= 3:
        return max(contigs[1:-1])
    else:
        return max(contigs)
score_dict['contig_not_hp_internal_max'] = score_dict['sequence'].map(seq_to_contig_not_hp_internal)

def seq_to_contig_not_hp_average(seq):
    nonhp = ''
    for aa in seq:
        if aa in 'FIWLVMY':
            nonhp += ' '
        else:
            nonhp += aa
    contigs = [len(x) for x in nonhp.split()]
    return np.average(contigs)
score_dict['contig_not_hp_avg'] = score_dict['sequence'].map(seq_to_contig_not_hp_average)


def seq_to_contig_not_hp_avg_norm(seq):
    nonhp = ''
    for aa in seq:
        if aa in 'FIWLVMY':
            nonhp += ' '
        else:
            nonhp += aa
    contigs = [len(x) for x in nonhp.split()]
    return np.average(contigs) / (float(len(seq)) / len([''] + [x for x in seq if x in 'FIWLVMY']))
score_dict['contig_not_hp_avg_norm'] = score_dict['sequence'].map(seq_to_contig_not_hp_avg_norm)



assert len(score_dict.index.values) == 1, len(score_dict.index.values)
score_dict['score_per_res'] = float(score_dict.total_score) / float(score_dict.n_res)
score_dict['fa_atr_per_res'] = float(score_dict.fa_atr) / float(score_dict.n_res)
score_dict['fa_rep_per_res'] = float(score_dict.fa_rep) / float(score_dict.n_res)
score_dict['net_atr_per_res'] = (float(score_dict.fa_atr) + float(score_dict.fa_rep)) / float(score_dict.n_res)
score_dict['net_sol_per_res'] = (float(score_dict.fa_sol) + float(score_dict.fa_elec)) / float(score_dict.n_res)
score_dict['net_atr_net_sol_per_res'] = (float(score_dict.net_atr_per_res) + float(score_dict.net_sol_per_res))


if 'buried_np' in score_dict: score_dict['buried_np_per_res'] = float(score_dict.buried_np) / float(score_dict.n_res)
if 'exposed_np' in score_dict: score_dict['exposed_np_per_res'] = float(score_dict.exposed_hydrophobics) / float(score_dict.n_res)
if 'buried_np' in score_dict and 'exposed_np' in score_dict: score_dict['buried_minus_exposed_per_res'] = (float(score_dict.buried_np) - float(score_dict.exposed_hydrophobics)) / float(score_dict.n_res)
#score_dict['exposed_polar_per_res'] = score_dict.exposed_polars / score_dict.n_res
#score_dict['sa_per_res'] = (score_dict.exposed_hydrophobics + score_dict.exposed_polars) / score_dict.n_res

#cut sites
def seq_to_tryp_cut_sites(sequence):
    tempsequence = sequence + 'X'
    tryp_cut_sites = 0
    for i in range(len(sequence)):
        if tempsequence[i] in 'RK' and tempsequence[i+1] != 'P': tryp_cut_sites += 1
    return tryp_cut_sites
score_dict['tryp_cut_sites'] = score_dict['sequence'].map(seq_to_tryp_cut_sites)

def seq_to_chymo_cut_sites(sequence):
    tempsequence = sequence + 'X'
    chymo_cut_sites = 0
    for i in range(len(sequence)):
        if tempsequence[i] in 'FYW' and tempsequence[i+1] != 'P': chymo_cut_sites += 1
    return chymo_cut_sites
score_dict['chymo_cut_sites'] = score_dict['sequence'].map(seq_to_chymo_cut_sites)

def seq_to_chymo_with_LM_cut_sites(sequence):
    tempsequence = sequence + 'X'
    chymo_cut_sites = 0.0
    for i in range(len(sequence)):
        if tempsequence[i] in 'FYW' and tempsequence[i+1] != 'P': chymo_cut_sites += 1.0
        if tempsequence[i] in 'LM' and tempsequence[i+1] != 'P': chymo_cut_sites += 0.5
    return chymo_cut_sites
score_dict['chymo_with_LM_cut_sites'] = score_dict['sequence'].map(seq_to_chymo_with_LM_cut_sites)

seq_to_cut_site_analysis = {}
def cut_site_analysis(sequence):
    tempsequence = sequence+'X'
    nearest_chymo_cut_to_Nterm = 1
    nearest_chymo_cut_to_Cterm = 1
    nearest_tryp_cut_to_Nterm = 1
    nearest_tryp_cut_to_Cterm = 1
    while (tempsequence[nearest_chymo_cut_to_Nterm-1] not in 'FYW' or tempsequence[nearest_chymo_cut_to_Nterm] == 'P')  and nearest_chymo_cut_to_Nterm < len(sequence):
        nearest_chymo_cut_to_Nterm += 1
    while (tempsequence[nearest_tryp_cut_to_Nterm-1] not in 'RK' or tempsequence[nearest_tryp_cut_to_Nterm] == 'P')  and nearest_tryp_cut_to_Nterm < len(sequence):
        nearest_tryp_cut_to_Nterm += 1
    while (tempsequence[-1 -nearest_chymo_cut_to_Cterm] not in 'FYW' or tempsequence[- nearest_chymo_cut_to_Cterm] == 'P')  and nearest_chymo_cut_to_Cterm < len(sequence):
        #print >> sys.stderr, tempsequence, nearest_tryp_cut_to_Nterm, nearest_chymo_cut_to_Nterm, nearest_chymo_cut_to_Cterm
        nearest_chymo_cut_to_Cterm += 1
    while (tempsequence[-1 -nearest_tryp_cut_to_Cterm] not in 'RK' or tempsequence[- nearest_tryp_cut_to_Cterm] == 'P')  and nearest_tryp_cut_to_Cterm < len(sequence):
        nearest_tryp_cut_to_Cterm += 1
    nearest_tryp_cut_to_term = min([nearest_tryp_cut_to_Cterm, nearest_tryp_cut_to_Nterm])
    nearest_chymo_cut_to_term = min([nearest_chymo_cut_to_Cterm, nearest_chymo_cut_to_Nterm])
    return (nearest_chymo_cut_to_Nterm, nearest_chymo_cut_to_Cterm, nearest_tryp_cut_to_Nterm, nearest_tryp_cut_to_Cterm, nearest_tryp_cut_to_term, nearest_chymo_cut_to_term)
for seq in score_dict.sequence.values:
    seq_to_cut_site_analysis[seq] = cut_site_analysis(seq)
score_dict['nearest_chymo_cut_to_Nterm'] = score_dict['sequence'].map(lambda x: seq_to_cut_site_analysis[x][0])
score_dict['nearest_chymo_cut_to_Cterm'] = score_dict['sequence'].map(lambda x: seq_to_cut_site_analysis[x][1])
score_dict['nearest_tryp_cut_to_Nterm'] = score_dict['sequence'].map(lambda x: seq_to_cut_site_analysis[x][2])
score_dict['nearest_tryp_cut_to_Cterm'] = score_dict['sequence'].map(lambda x: seq_to_cut_site_analysis[x][3])
score_dict['nearest_tryp_cut_to_term'] = score_dict['sequence'].map(lambda x: seq_to_cut_site_analysis[x][4])
score_dict['nearest_chymo_cut_to_term'] = score_dict['sequence'].map(lambda x: seq_to_cut_site_analysis[x][5])

seq_column = score_dict.sequence
desc_column = score_dict.description
score_dict = score_dict.drop('sequence', 1)
score_dict = score_dict.drop('description', 1)
score_dict['sequence'] = seq_column
score_dict['description'] = desc_column

# Output the results in `json` format as a dictionary, outputting it to the
# terminal using sys.stdout
print >> sys.stderr, "Concluded enhance_scorefile.py and will return results"
json.dump(fp=sys.stdout, obj=score_dict.to_dict(orient="records"), indent=2)
