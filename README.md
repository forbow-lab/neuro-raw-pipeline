# neuro-raw-pipeline
Raw data ingestion and organization

### Download DICOMS (T1w/T2w/DWI)

- Prior to leaving 3T MRI request the MRI Technologist "push" the exam DICOMS to the PBIL-master AE Node. 
- From an iMac in the PBIL, download the DICOMS from PBIL-master onto the Shared drive:
  - Open Finder window, create a new subject folder on shared drive (e.g. `/shared/uher/FORBOW/rawdata/101_C_20190225/`)
    - Open Horos
    - Select "Query" on main interface
        - Search "Forbow", select date-range (e.g. "last 7 days")
        - Select "Query" from middle of query window, within seconds a list of all Forbow scans in last 7 days should appear.
        - Download new exams by hitting the green download button beside each FORBOW participant
        - Monitor the "Activity" panel on the left side of the main Horos window until all downloads are complete.
        - Before closing the query box, ensure that no files are still transferring in the "Activity" panel. 
        - From main Horos window, select four main sequences [T1w, T2w, DTI-30DIR, DTI_B0] for one participant, and drag into Finder window to drop on the subject folder.
        - Monitor the "Activity" panel again until all DICOMS are exported to the subject folder on the shared drive.
        - Rename and organize the exported DICOM folder under the subject folder (e.g. `/shared/uher/FORBOW/rawdata/101_C_20190225/DICOMS/`).
        - Note, if there are multiple T1w_BRAVO and/or T2w_FlairPrep sequences, only copy over the highest quality scan.

---

### Download Pfiles (RS)

- Confirm Pfiles are ready to download from remote BIOTIC server: 
* `ssh biotic@lauterbur` 
  * enter password
* Navigate to our data folder on their server `cd /biotic/data/Forbow/`
* `ls` to see the list of participants
* `ls` into specific participant folder
* If participant files appear with a timestamp they are good to go, if not, you'll have to wait for BIOTIC to convert them 
* `exit` to close connection


After the EPI images have been transferred to the server at BIOTIC (lauterbur):
* open terminal and `cd` into SUBJECT ID in `rawdata folder`
* run the following command *all one line* replacing IDs with the correct subject IDs and date:
* `nohup rsync -a -e "ssh" biotic@lauterbur:/biotic/3Tdata/Forbow/###_X/ /shared/uher/FORBOW_Brain/neuro_data/Biotic3T/rawdata/###_X_YYYYMMDD/RS/`
  * copies the resting state data from the lauterbur server to our server into the participants RS rawdata folder
  * make sure to supply correct subject ID and scan label for `###_X` and the correct date of scan for `YYYYMMDD`
  * enter password

Note: may take 10-15 minutes to copy across network (~17 GB)


---

### Conversion: Resting State
(typically 1-3 hours per subject)
The resting state sequence from the scanner is custom built for us as of now, and as such requires a number of extra steps. One of which was done above when SSH transferring the EPI sequences from lauterbur to PBIL separately. Once the data is full transferred, we need to transfer the raw p files into a format we can use: nifti (.nii)

There should be 3 sets of p files, this is how they correspond to the scanner sequences:


Scanner     			    | 	Time   | Pfile | Description
--- 	 				      	| 	---    | ---    | ---
EPI NoMUX Cal		     	|	00:28    | ~700 mb, lowest numeric title    |   This is a scan with no multiband factor and minimal distortion, will later be used as an SBRef in the minimal preprocessing pipelines.
EPI Bomap - Rev		    |	00:15    | ~300mb, higher numeric title   | This is a few volumes of the RS sequence acquired in a reversed phase encode direction. Will be used to estimate and correct the distortion.
EPI MUX 3MM RS		    |	08:06   | ~15gb, highest numeric title    | This is the main resting state scan where we acquire 500 volumes in a sub-second TR. As you can see it is by far the largest, and if it is smaller than 15gb then something went wrong with the scan.


Note, each Pfile is also accompanied by `_noise` `_param` and `_ref` dat files. There may be a lone Pfile, lowest in numeric title and a small size, unaccompanied by the above dat files. This is likely the shim, so **create a new folder** in the subject RS directory called `shim` and **move** that file there.
  * Update: now this is automatically taken care of in the conversion script. The p file could also be spectroscopy data from other studies caught by a wildcard `*` so for example `P*` catches anything that starts with a P and ends with whatever the rest of the string might be.


Now that the `rawdata` subject RS folder is organized, we will convert the RAW p files into nifti, and compress the original raws.

1. Open the terminal.
2. Drag the `run_mux2niismarter.sh` script into the terminal.
3. Then drag the subject `RS` folder into the same line in the terminal.
4. Hit enter and wait for the conversion to complete before moving on to the step.

Note this might take some time, and you can do multiple subjects sequentially by dragging additional RS folders into the line with the `run_mux2nii.sh` script. Likewise you can do multiple subjects in parallel by opening new terminal windows or using the GNU parallel code found in the LGI pipeline.

---

### Conversion: Everything else
(typically ~10 minutes per subject)

This step assumes that you have finished: transferring T1w, T2w, DTI data from OsiriX, SSH transferring RS data, and converting RS data to nifti.

The following steps convert the dicom data in the `rawdata` folder into nifti format and copy the nifti data up one directory into `Biotic3T` and create and organize the subject folders in a specific way to make later preprocessing possible.

1. Open the terminal.
2. Navigate to the `_scripts` folder in /shared/uher/FORBOW_Brain/neuro_data/Biotic3T/
3. Drag the `0_convert_dcm2niix.sh` script into the terminal
  * alternatively, cd into the `_scripts` folder and type `0` and click **tab** to autocomplete
4. In terminal, after the script path, type subject IDs that need to be converted (space separated), for example `032_A 033_A 031_B`
5. Click **enter** to run script on specified subjects.

Check /shared/uher/FORBOW_Brain/neuro_data/Biotic3T/ for the subjects specified. There should be new folders with the subject IDs. In those folders you will find the `unprocessed` folder which has all the nifti files ready for preprocessing, and a `log` folder where preprocessing pipeline will log all the steps and spit out any errors.

---
