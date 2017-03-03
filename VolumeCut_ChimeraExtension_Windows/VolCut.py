import chimera
from chimera import runCommand
import sys
import subprocess
import os

# Cuts the volume and then opens up the new cut to visualize it. 
# TODO: Add more error handling. Most is done on the gui side. 
def cutvol(pdb_loaded_num, mrc_loaded_num, radius, resolution, chain):
    # get the path of both the MRC and PDB files. 
    path_pdb = str(chimera.openModels.list()[pdb_loaded_num].openedAs[0])
    path_mrc = str(chimera.openModels.list()[mrc_loaded_num].openedAs[0])

 # # Define the path to the executable.
# path_to_volcut = str(os.path.realpath(__file__))
# path_to_volcut = path_to_volcut[0:len(path_to_volcut) - 3]
# print path_to_volcut
# path_to_volcut = "\"" + path_to_volcut + ".exe\""
# print path_to_volcut
# type(path_to_volcut)
# # Set up args to pass to binary
# args = (str(path_to_volcut), path_pdb, path_mrc, str(radius), str(resolution), chain)

	# Run the binary to cut file, wait for it to finish, print output
    popen = subprocess.Popen(args, stdout=subprocess.PIPE)
	
    # Set up args to pass to binary
    args = ("VolCut.exe", path_pdb, path_mrc, str(radius), str(resolution), chain)

	# Run the binary to cut file, wait for it to finish, print output
    popen = subprocess.Popen(args, stdout=subprocess.PIPE)
    popen.wait()
    output = popen.stdout.read()
    print 'cutting with parameters: ' + output
	
    # Get the new file names as they're written in the binary. 
    cut_mrc_name = chimera.openModels.list()[mrc_loaded_num].name[0:len(chimera.openModels.list()[mrc_loaded_num].name) - 4] + '_' + chain + '_cut.mrc'
    cut_pdb_name = chimera.openModels.list()[pdb_loaded_num].name[0:len(chimera.openModels.list()[pdb_loaded_num].name) - 4] + '_' + chain + '_cut.pdb'
    
	# Close any open models if the name is the same (otherwise chimera will just open the same file)
    for model in chimera.openModels.list():
        if cut_mrc_name == model.name or cut_pdb_name == model.name:
            runCommand('close' + str(model))
	
	# Open the cut PDB and MRC files so the user can save them if they want
    opencmdmrc = 'open ' + os.getcwd() + '/' + cut_mrc_name
    opencmdpdb = 'open ' + os.getcwd() + '/' + cut_pdb_name
    runCommand(opencmdmrc)
    runCommand(opencmdpdb)