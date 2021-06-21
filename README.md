# neuro-raw-pipeline
Raw data ingestion and organization

### 1. Download DICOMS (T1w/T2w/DWI)

- Prior to leaving 3T MRI request the MRI Technologist "push" the exam DICOMS to the PBIL-master Node. 
- From an iMac in the PBIL, download the DICOMS from PBIL-master onto the Shared drive:
  * Open Finder window, create a new subject folder on shared drive (e.g. `/shared/uher/FORBOW/rawdata/100_C_20190320/`)
  * Open Horos
    - Select "Query" on main interface
    - Search "Forbow", select date-range (e.g. "last 7 days")
    - Select "Query" from middle of query window, within seconds a list of all Forbow scans in last 7 days should appear.
    - Download new exams by hitting the green download button beside each FORBOW participant
    - Monitor the "Activity" panel on the left side of the main Horos window until all downloads are complete.
    - Before closing the query box, ensure that no files are still transferring in the "Activity" panel. 
    - From main Horos window, select four main sequences [T1w, T2w, DTI-30DIR, DTI_B0], drag into Finder window to drop on subject folder.
    - Monitor the "Activity" panel again until all DICOMS are exported to subject folder on shared drive.
    - Rename and organize the exported DICOM folder under subject folder (e.g. `/shared/uher/FORBOW/rawdata/100_C_20190320/DICOMS/`).
    - Note, if there are multiple T1w_BRAVO and/or T2w_FlairPrep sequences, only copy over the highest quality scan.

---

### 2. Data Quality ratings and Reminders

- Create a  "readme" file that contains pertinent notes about that particular participant and scan (e.g. T1 BRAVO collected twice due to motion, kept 2nd ) as well as the exam number
- Create a file to rate the quality of the T1 BRAVO, T2 CUBE, and T2 PROMO on a four point scale (e.g. 109_A_qc_ratings). See any previous rawdata folder for an example
- Complete the remaining fields of the neurolog
- Send images to participants that requested a picture of their brain 

---

### 3. Download Pfiles (RS)

- First. confirm Pfiles are ready to download from remote BIOTIC server: 
* `ssh biotic@lauterbur` (enter password)
* Navigate to data folder on remote server `cd /biotic/data/Forbow/`
* `ls` to see the list of participant folders
* `ls 100_C` to list contents of specific participant folder
* If participant files appear with timestamps they are ready to download, if not, you'll have to wait for BIOTIC to convert them. 
* `exit` to close remote connection.

- Once Pfiles are ready, download from remote BIOTIC server:
* open terminal and `cd` into SUBJECT ID (e.g., `/shared/uher/FORBOW/rawdata/100_C_20190320/`
* use rsync command below (*all one line*) replacing IDs with the correct subject IDs and date:
* `rsync -av biotic@lauterbur:/biotic/data/Forbow/100_C/ /shared/uher/FORBOW/rawdata/100_C_20190320/RS/`
  * this copies RS Pfiles from remote server to our local server into the participant's RS folder.
  * make sure to supply correct subject ID and scan label for `###_X` and the correct date of scan for `YYYYMMDD`
  * enter password if requested, watch terminal to monitor file transfer until complete and returned to prompt.

Note: may take 5-10 minutes to copy across network (~17 GB)


---


### 4. Run Rawdata Pipeline Script to convert, deface, and QC.

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
  * start pipeline script: `nohup ./_1_scripts/FORBOW_0_Run_RAWDATA_Pipeline.sh 101_C  >nohup_rawpipe_100_C_20190411.txt 2>&1 &`
  * monitor script is running with ps command: `ps aux | grep '100_C'` or htop command ('q' to quit htop)
  * pipeline takes 2-3 hrs to complete (automatically uses 4-8 CPU Cores on Jaylah)
  * from Finder window, check subject _QC_report.html to confirm completeness, eg: `open /shared/uher/FORBOW/rawdata/100_C_20190320/NIFTIS/QC/_QC_report.html`


---
Note, the Rawdata Pipeline script runs each of the below scripts, in order:
  1. Convert raw Pfiles [1_run_mux2nii.sh](https://github.com/forbow-lab/neuro-raw-pipeline/wiki/GE-3T-EPI-Pfile-Conversion-\(mux2nii\))
  2. Convert DICOMS into NIFTIS [2_create_niftis.sh](https://github.com/forbow-lab/neuro-raw-pipeline/wiki/DICOM-to-NIFTI-Conversion)
