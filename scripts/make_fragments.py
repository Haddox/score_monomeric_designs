"""
Python script for making and analyzing fragments for an input protein design
"""


from argparse import ArgumentParser
import os
import glob
from multiprocessing import Pool
parser  = ArgumentParser()
import sys
from shutil import copy

fragment_tools = '/work/robetta/workspace/labFragPicker_DO_NOT_REMOVE/Rosetta/tools/fragment_tools/'
psipred = '/work/robetta/workspace/labFragPicker_DO_NOT_REMOVE/psipred3.21/'
scripts_dir = '/work/robetta/workspace/labFragPicker_DO_NOT_REMOVE/bakerlab_scripts/boinc/'
nnmake = '/work/robetta/workspace/labFragPicker_DO_NOT_REMOVE/nnmake/pNNMAKE.gnu'
csbuild = '/work/robetta/workspace/labFragPicker_DO_NOT_REMOVE/csbuild/'
cm_scripts = '/work/robetta/workspace/labFragPicker_DO_NOT_REMOVE/cm_scripts/bin/'
rosetta = '/work/robetta/workspace/labFragPicker_DO_NOT_REMOVE/Rosetta/main/source/bin/'
if not 'SPARKSDIR' in os.environ:
  os.environ[ 'SPARKSXDIR' ] = '/work/robetta/workspace/labFragPicker_DO_NOT_REMOVE/sparks-x'

def chckmkdir(path,out = sys.stdout):
  if not os.path.exists(path):
    os.mkdir(path)
    return True
  else:
    return False

def generate_fragments(pdb, args):
  currentdir = os.getcwd()
  basename, extension = os.path.splitext(pdb)
  pdbname = os.path.split(basename)[-1]
  exec_dir = currentdir +"/" +  pdbname + "_fragments"
  exec_dirs = []
  dir_maked  = chckmkdir(exec_dir)
  if dir_maked:
    copy(pdb, exec_dir + "/00001.pdb")
    os.chdir(exec_dir)
    os.system(cm_scripts + "pdb2fasta.py  00001.pdb >  00001.fasta")
    if (args.single_psipred and (not args.nohoms)):
      if args.csbuild:
        os.system(psipred + "runpsipred_csbuild_single 00001.fasta")
        os.system(csbuild+"csbuild -i  00001.fasta -I fas -D "+csbuild+"csblast-2.2.3/data/K4000.crf -o  00001.check -O chk")
        if args.v2012_db:
          os.system(fragment_tools + "make_fragments_tjTest.pl  -nosam -id 00001 00001.fasta -psipredfile 00001.csb.ss2 -nopsiblast_profile")
        else:
          os.system("cp 00001.csb.ss2 00001.psipred_ss2")
          copy(scripts_dir +"/path_defs.txt", exec_dir)
          os.system(scripts_dir+"/make_checkpoint.pl -id 00001 -fasta 00001.fasta")
          os.system(nnmake + " aa 0000 1")
          os.system("mv aa0000103_05.200_v1_3 00001.200.3mers")
          os.system("mv aa0000109_05.200_v1_3 00001.200.9mers")
      else:
        os.system(psipred + "runpsipred_single 00001.fasta")
        os.system(fragment_tools + "make_fragments_tjTest.pl  -nosam -id 00001 00001.fasta -psipredfile 00001.ss2 -nopsiblast_profile")
    else:
      os.system(fragment_tools + "make_fragments.pl  -nosam -nohoms -id 00001 00001.fasta")
    if(args.eval_method == 'nobu'):
      os.system(rosetta+"r_frag_quality.linuxgccrelease -in:file:native 00001.pdb -f  00001.200.9mers")
      os.system(scripts_dir+"/count_goodfrag.pl --s=frag_qual.dat > good_fragments_count.txt")
      os.system('python2 /work/grocklin/gabe_scripts/frag_qual_and_sequence.py frag_qual.dat')
    if(args.eval_method == 'tj'):
      os.system(rosetta+"fragment_rmsd.linuxgccrelease -in:file:s 00001.pdb -in:file:frag3 00001.200.9mers > rmsd_9mers.txt")
  os.chdir(currentdir)

def generate_fragments_single_arg(pdb_args):
    generate_fragments(pdb_args[0],pdb_args[1])
#Note: the multiprocessor support for pool takes only a single argument in python 2.7. If we were running python 3.3 I believe starmap could be used making this function unncessary


if __name__ == '__main__':
  parser.add_argument('-n_proc', default=10, type=int, help='number of cores to use, by default 10')
  parser.add_argument('-single_psipred', default = True,  help='use psipred single sequence prediction')
  parser.add_argument('-no-single_psipred',dest='single_psipred',action='store_false')
  parser.add_argument('-nohoms', default = False,action='store_true', help='uses BLAST to detect homologs')
  parser.add_argument('-csbuild',default= True ,help='use csbuild to generate the checkpoint file only valid for single sequence')
  parser.add_argument('-no-csbuild',dest='csbuild',action='store_false')
  parser.add_argument('-pdbs', nargs='+', help='pdbs to proccess')
  parser.add_argument('-pdbs_folder', type=str, help=' directory where the pdbs to process exist. (Must be labeled *.pdb) ')
  parser.add_argument('-eval_method',choices=['nobu','tj','none'], default='nobu', help='evaluation method <nobu>,<tj>,<none> only valid options')
  parser.add_argument('-v2012_db', default= False, action='store_true', help='uses 2012 + C++ frag picker')

  args = parser.parse_args()

  pdbs = glob.glob(str(args.pdbs_folder)+"/*.pdb")[:2]
  if (len(pdbs) == 0):
    if(args.pdbs):
      pdbs = args.pdbs
    else:
      print("NO PDBS ENTERED")
      sys.exit()

  cmd_list = [( [pdb,args] ) for pdb in pdbs ]
  myPool    = Pool( processes = args.n_proc )
  myResults = myPool.map( generate_fragments_single_arg,cmd_list )
