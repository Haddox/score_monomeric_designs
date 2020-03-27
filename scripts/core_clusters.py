#!/usr/local/bin/python
import sys
import Bio.PDB
import networkx as nx
from numpy import *
from Bio.PDB import *

def main( files ):
    hydrophobic_names=['PHE','ILE','LEU','MET','VAL','TRP','TYR']#,'ALA']

    canon="""ALA
ARG
ASN
ASP
CYS
GLN
GLU
GLY
HIS
ILE
LEU
LYS
MET
PHE
SER
THR
TRP
TYR
VAL
PRO""".split('\n')
    output = ""
    for (filenum, fn) in enumerate(files):
        structure = Bio.PDB.PDBParser(QUIET=True).get_structure(fn,fn)
        contacts={}
        model = structure[0]
        polypeptides = Bio.PDB.PPBuilder().build_peptides(model)
        seq = ''
        for p in polypeptides:
            seq += p.get_sequence()


        residues = [res for res in model.get_residues()]
        all_canonical = [res.get_resname() in canon for res in residues]
        G = nx.Graph()
        if False in all_canonical: continue

        for i, res in enumerate(residues):
            if res.get_resname() in hydrophobic_names:
                G.add_node(i+1)

        for i, res in enumerate(residues):
            for j, res2 in enumerate(residues):
                if res.get_resname() not in hydrophobic_names or res2.get_resname() not in hydrophobic_names: continue
                if i + 2 < j:
                    contact = False
                    for a1 in [x for x in res if ('H' not in x.get_name() and x.get_name() not in ['C','CA','O','N'])]:
                        for a2 in [x for x in res2 if ('H' not in x.get_name() and x.get_name() not in ['C','CA','O','N'])]:
                            if a1 - a2 < 4.3:
                                G.add_edge(i+1,j+1)
                                #print i+1, j+1

                                continue

        #print G.nodes()
        #print G.edges()
        subgraphs = list( nx.connected_components(G) )

        output += fn + ' ' + str(seq) + ' '
        output += str( max([0] + [len(x) for x in subgraphs]) ) + ' '
        output += "%d " % len(subgraphs)
        output += "%d " % len(G.edges())
        output += "%.4f " % ( float(len(G.edges()))/max([1,float(len(G.nodes()))]) )
        output += subgraphs_to_str( subgraphs ) + '\n'
        if filenum > 0 and filenum % 100 == 0:
            print >>sys.stderr, "Processed", filenum, "/", len(files), "structures."
    return output.rstrip( '\n' )

def subgraphs_to_str(x):
    out = ''
    for subgraph in x:
        out += '+'.join('%s' % x for x in sorted(subgraph))
        out += '_'
    return out[0:-1]

if __name__ == "__main__":
    print main( sys.argv[1:] )
