# neuro-raw-pipeline
Raw data ingestion and organization

### Download DICOMS (T1w/T2w/DWI)

- Prior to leaving 3T MRI request the MRI Technologist "push" the exam DICOMS to the PBIL-master Node. 
- From an iMac in the PBIL, download the DICOMS from PBIL-master onto the Shared drive:
  * Open Finder window, create a new subject folder on shared drive (e.g. `/shared/uher/FORBOW/rawdata/101_C_20190225/`)
  * Open Horos
    - Select "Query" on main interface
    - Search "Forbow", select date-range (e.g. "last 7 days")
    - Select "Query" from middle of query window, within seconds a list of all Forbow scans in last 7 days should appear.
    - Download new exams by hitting the green download button beside each FORBOW participant
    - Monitor the "Activity" panel on the left side of the main Horos window until all downloads are complete.
    - Before closing the query box, ensure that no files are still transferring in the "Activity" panel. 
    - From main Horos window, select four main sequences [T1w, T2w, DTI-30DIR, DTI_B0], drag into Finder window to drop on subject folder.
    - Monitor the "Activity" panel again until all DICOMS are exported to subject folder on shared drive.
    - Rename and organize the exported DICOM folder under subject folder (e.g. `/shared/uher/FORBOW/rawdata/101_C_20190225/DICOMS/`).
    - Note, if there are multiple T1w_BRAVO and/or T2w_FlairPrep sequences, only copy over the highest quality scan.

---

### Download Pfiles (RS)

- First. confirm Pfiles are ready to download from remote BIOTIC server: 
* `ssh biotic@lauterbur` (enter password)
* Navigate to data folder on remote server `cd /biotic/data/Forbow/`
* `ls` to see the list of participant folders
* `ls 101_C` to list contents of specific participant folder
* If participant files appear with timestamps they are ready to download, if not, you'll have to wait for BIOTIC to convert them. 
* `exit` to close remote connection.


- Once Pfiles are ready, download from remote BIOTIC server:
* open terminal and `cd` into SUBJECT ID (e.g., `/shared/uher/FORBOW/rawdata/101_C_20190225/`
* use rsync command below (*all one line*) replacing IDs with the correct subject IDs and date:
* `rsync -av biotic@lauterbur:/biotic/data/Forbow/###_X/ /shared/uher/FORBOW/rawdata/###_X_YYYYMMDD/RS/`
  * this copies RS Pfiles from remote server to our local server into the participant's RS folder.
  * make sure to supply correct subject ID and scan label for `###_X` and the correct date of scan for `YYYYMMDD`
  * enter password if requested, watch terminal to monitor file transfer until complete and returned to prompt.

Note: may take 5-10 minutes to copy across network (~17 GB)


---


### Run Rawdata Pipeline Script to convert, deface, and QC.

- The Rawdata Pipeline Script:
  1. converts RS Pfiles into NIFTI format (1-2 hrs per subject),
  2. DICOM sequences into NIFTI,
  3. structurally defaces the anatomicals,
  4. generates a simple QC HTML page for ensuring sequences are ready for analysis, 
  5. creates BIDS folder for functional processing.

- From McCoy/Picard iMac:
  * Open Terminal
  * remote login to Jaylah: `ssh uher@jaylah`
  * navigate to rawdata folder: `cd /shared/uher/FORBOW/rawdata/`
  * start pipeline script: `nohup ./_1_scripts/FORBOW_0_Run_RAWDATA_Pipeline.sh 101_C  >nohup_rawpipe_101_C_20190312.txt 2>&1 &`
  * monitor script is running with ps command: `ps aux | grep '101_C'` or htop command ('q' to quit htop)
  * pipeline takes 2-3 hrs to complete (automatically uses 4-8 CPU Cores on Jaylah)
  * from Finder window, check subject _QC_report.html to confirm completeness, eg: `open /shared/uher/FORBOW/rawdata/101_C_20190225/NIFTIS/QC/_QC_report.html`


---

### RS Pfile to NIFTI Conversion Details:
(typically 1-3 hours per subject)
The resting state sequence from the scanner is custom built for us, and as such requires a number of extra steps. Manually downloading these is one of the extra steps. Once the data is fully transferred, a second step is requied to 'recon' the raw Pfiles (conversion from frequency "K-space" with 32-channels and 3x multi-banding into image-space), and then saved into neurimaging data format NIFTI (.nii/.nii.gz).

There should be 3 sets of Pfiles, this is how they correspond to the scanner sequences:


Scanner     			    | 	Time   | Pfile | Description
--- 	 				      	| 	---    | ---    | ---
EPI NoMUX Cal		     	|	00:28    | ~700 mb, lowest numeric title    |   This is a scan with no multiband factor and minimal distortion, will later be used as an SBRef in the minimal preprocessing pipelines.
EPI Bomap - Rev		    |	00:15    | ~300mb, higher numeric title   | This is a few volumes of the RS sequence acquired in a reversed phase encode direction. Will be used to estimate and correct the distortion.
EPI MUX 3MM RS		    |	08:06   | ~15gb, highest numeric title    | This is the main resting state scan where we acquire 500 volumes in a sub-second TR. As you can see it is by far the largest, and if it is smaller than 15gb then something went wrong with the scan.


Note, each Pfile is also accompanied by `_noise.dat` `_param.dat` and `_ref.h5` files. There may be a lone Pfile, lowest in numeric title and a small size, unaccompanied by the above dat files. This is likely the shim, so **create a new folder** in the subject RS directory called `shim` and **move** that file there.
  * Update: now this is automatically taken care of in the conversion script. The p file could also be spectroscopy data from other studies caught by a wildcard `*` so for example `P*` catches anything that starts with a P and ends with whatever the rest of the string might be.


Now that the `rawdata` subject RS folder is organized, we will convert the RAW p files into nifti, and compress the original raws.

1. Open the terminal.
2. Drag the `run_mux2niismarter.sh` script into the terminal.
3. Then drag the subject `RS` folder into the same line in the terminal.
4. Hit enter and wait for the conversion to complete before moving on to the step.

Note this might take some time, and you can do multiple subjects sequentially by dragging additional RS folders into the line with the `run_mux2nii.sh` script. Likewise you can do multiple subjects in parallel by opening new terminal windows or using the GNU parallel code found in the LGI pipeline.

---

### DICOM to NIFTI Conversion Details:
(typically ~10 minutes per subject)

This step assumes the T1w, T2w, DTI data are downloaded to the shared drive from Horos, and RS Pfiles have been downloaded, organized, and converted into NIFTI files.

The following steps convert the dicom data in the `rawdata` folder into nifti format and copy the nifti data up one directory into `Biotic3T` and create and organize the subject folders in a specific way to make later preprocessing possible.

1. Open the terminal.
2. Navigate to the `_scripts` folder in /shared/uher/FORBOW_Brain/neuro_data/Biotic3T/
3. Drag the `0_convert_dcm2niix.sh` script into the terminal
  * alternatively, cd into the `_scripts` folder and type `0` and click **tab** to autocomplete
4. In terminal, after the script path, type subject IDs that need to be converted (space separated), for example `032_A 033_A 031_B`
5. Click **enter** to run script on specified subjects.

Check /shared/uher/FORBOW_Brain/neuro_data/Biotic3T/ for the subjects specified. There should be new folders with the subject IDs. In those folders you will find the `unprocessed` folder which has all the nifti files ready for preprocessing, and a `log` folder where preprocessing pipeline will log all the steps and spit out any errors.

---
