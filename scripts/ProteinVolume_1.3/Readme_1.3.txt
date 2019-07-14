==================================================================================================================
ProteinVolume 1.3 Readme
By Calvin Chen, Rensselaer Polytechnic Institute (RPI) - chenc13@rpi.edu
Last Updated 1/10/2014
==================================================================================================================
Requirements: Java Runtime Environment 1.7 and above

Usage: java -jar ProteinVolume.jar [options] filename/directory

Flags:
  -het,  --heteroatoms [off]                      Read heteroatoms
  -rf,   --radiusFileName [bondi.rad]             Use specified van der Waals radius file. Setting to a nonexistent file will result in a very sparse hard-coded radius set.
  -lf,   --listFileName [off]                     Use specified protein list file

  -sr,   --surfaceResolution [64]                 Set surface angle increment between surface probes (PI/#)
  -vpr,  --volumeProbeRadius [0.08 0.02]          Set starting and smallest volume probe radii (angstrom)
  -wpr,  --waterProbeRadius [1.4]                 Set water radius used to generate molecular surface (angstrom)
  -spmd, --surfaceProbeMinDistance [0.05]         Set minimum distance between surface probes (angstrom)

  -v,    --version,                               Print version
  -d,    --debug,                                 Print time taken to calculate individual components and # of probes
  
==================================================================================================================
How to run from command line:

[path to java] [java flags] -jar ProteinVolume.jar [ProteinVolume flags] "[list of protein directories and/or filenames [||starting protein #||ending protein #] ]"

Note: Surround all paths and range numbers with double quotes. "/data/test/1UBQ||100||200" (Proteins 100-200 in /data/test/1UBQ directory)

A template run script has been included (runProteinVolume.bat/.sh)

Running on manually specified protein directories or protein filenames:
1. Edit runProteinVolume.bat/.sh
2. List protein directories and/or filenames you wish to analyze (e.g "java -jar C:\ProteinDir C:\1UBQ.pdb C:\2ACY.pdb". All .pdb files in C:\ProteinDir will be processed along with C:\1UBQ.pdb and C:\2ACY.pdb)
4. Run runProteinVolume.bat/.sh

Running on proteins specified in a list file:
1. Open your list file (text format)
2. Enter protein directories and/or filenames you wish to analyze separated with line break (Enter)
3. Add "-f" flag to the ProteinVolume command (e.g. "java -jar ProteinVolume.jar -f")
4. Specify path of your list file (e.g "java -jar ProteinVolume.jar -f C:\listFile.txt")
5. Run runProteinVolume.bat/.sh

To specify numerical range of proteins to process in a protein directory:
Specify the starting and ending protein (optional) number after the protein directory (e.g. "C:\ProteinDir||20||50". In this case, it will only process proteins 20-50 in that directory). If the starting protein number is 0, it will process starting from the first protein in the directory. If the ending protein number is not specified or is 0, it will process up til the last protein in that directory.
==================================================================================================================
Custom radii:
You may specify any radius value for any element type by creating a text file that is space or tab delimited e.g.:
C	1.7
N	1.55
H	1.2
and so forth

Setting H radius to 0 or omitting H entirely from the radius file will tell PV to ignore hydrogens from the PDB file.
==================================================================================================================
Where output is stored:

1. Console output
2. Output (%date%).txt (in the specified protein directory)
==================================================================================================================
Output format:

Version: [version]
Date: [date]
Protein directory: [protein directory]
Volume probe radius: [starting probe size] -> [ending probe size]
Surface probe minimum distance: [surface probe minimum distance]
Starting and ending indices: [start index] to [end index] (if indices are specified)
Columns: Protein Name, Total Volume, Void Volume, VDW Volume, Packing Density, Time Taken (ms)

Sample output:
PV Version: 1.3
Date: 2017-01-10 02:45:25 PM
Protein: [Path]\1UBQ_001.pdb
Volume probe radius: 0.080A -> 0.020A
Surface probe minimum distance: 0.1A
Number of atoms and residues: 1231, 76

Protein               Total Volume (A3)      Void Volume       VDW Volume  Packing Density  Time Taken (ms)
1UBQ_001                      10230.378         2483.516         7746.862            0.757  27105