# Volume Cut and Volume Compare
Currently contains the source code for Volume Cut and Volume Compare. 
Compiled Source Code can be found in /LinuxCompiled folder.

UCSF Chimera extension module can be found in the VolumeCut_ChimeraExtension_Windows folder. Below are step-by-step installation instructions:

1) 
Create a new folder, then place the VolumeCut_ChimeraExtension_Windows folder from the repository into it.

2)
Open UCSF Chimera and select Favorites -> Add to Favorites/Toolbar
![alt tag](https://github.com/nissensonm/VolumeCut/blob/master/img/step1.jpg)

3) 
Click the add button in the tools menu.
![alt tag](https://github.com/nissensonm/VolumeCut/blob/master/img/step2.jpg)

4)
In the "Add Extension Directory" menu, navigate to the folder one level above the VolumeCut_ChimeraExtension_Windows folder. Do not select the VolumeCut_ChimeraExtension_Windows folder itself, rather you must select the folder in the file structure one level above.
![alt tag](https://github.com/nissensonm/VolumeCut/blob/master/img/step3.jpg)

5)
After selecting the folder one level above the VolumeCut_ChimeraExtension_Windows folder, select save.
![alt tag](https://github.com/nissensonm/VolumeCut/blob/master/img/step4.jpg)

6) 
Back on the UCSF Chimera menu, navigate to Tools, select Utilities, then select VolumeCut.
![alt tag](https://github.com/nissensonm/VolumeCut/blob/master/img/step5.jpg)

7)
The volume cut window should now appear. They will be blank until a MRC and PDB are loaded. 
![alt tag](https://github.com/nissensonm/VolumeCut/blob/master/img/step6.jpg)

8)
Load a MRC file and fitted PDB file in UCSF Chimera. For this example, I used EMDB-2513 and PDB-4ci0.
Fill out the Resolution of the MRC file (3.36 is the stated resolution) and desired Radius to cut. A radius of 2.5 was selected as it selects most of the information from a chain without cutting too much voxel data.
![alt tag](https://github.com/nissensonm/VolumeCut/blob/master/img/step7.jpg)

9)
Once all information has been filled out and selected, click the Cut Volume button. This then isolates the PDB chain selected, the voxels surrounding the chain, and then loads the results in UCSF Chimera. Below the original volume voxel data has been made transparent to highlight the cut volume data (in yellow).
![alt tag](https://github.com/nissensonm/VolumeCut/blob/master/img/step8.jpg)

10) 
To save the cut MRC, select the cut volume in the volume viewer, then go to File -> Save map as...
![alt tag](https://github.com/nissensonm/VolumeCut/blob/master/img/step9.jpg)

11) 
To save the cut PDB, navigate to the main UCSF Chimera tool-bar menu, navigate to File -> Save PDB...
![alt tag](https://github.com/nissensonm/VolumeCut/blob/master/img/step10.jpg)
