function d = mux2nii(matfile)
%% mux2nii
%
% Convert processed mux epi data .mat file into .nii file
%
% mux2nii(matfile);
%
% Input:
%	matfile : 	Name of .mat file from mux epi processing (with .mat)
%
% Output:
%	None explicitly. Creates a Nifti file called {pfile}.nii (ie.
%	P12345.7.nii)

[pfilePath pfileName pfileExt] = fileparts(matfile);

% Read in .mat file
  if(exist(matfile, 'file') == 2)
    load(matfile)
  else
    error('Matlab file %s with processed epi mux data not found.',matfile);
  end

% Read in pfile header using GE Orchestra package
  pfileHandle = GERecon('Pfile.Load', pfile);
  header = GERecon('Pfile.Header', pfileHandle);

% Write nii file
  d = d(:,:,:,3:end);   % Discard first 2 reference volumes (possibly only needed for internal ref scans)
  d = flipdim(d,2);     % Reorient matrix so correct nii orientation
  nii = make_nii(d);

% Amend template header to account for acquired scan
  nii.hdr.dime.pixdim(2)    = header.ImageData.pixsize_X;   % X resolution
  nii.hdr.dime.pixdim(3)    = header.ImageData.pixsize_Y;   % Y resolution
  nii.hdr.dime.pixdim(4)    = header.ImageData.slthick;     % Z resolution
  nii.hdr.dime.pixdim(5)    = header.ImageData.tr/1000000;  % TR resolution
/* >> header.ImageData 
	...
     dim_X: 72
	 dim_Y: 72
 pixsize_X: 3
 pixsize_Y: 3
	 ctr_R: -3.9513
	 ctr_A: 32.9297
	 ctr_S: -41.5775
	norm_R: 0.0188
	norm_A: -0.1687
	norm_S: 0.9855
	tlhc_R: 107.8988
	tlhc_A: 135.7799
	tlhc_S: -26.1049
	trhc_R: -107.9191
	trhc_A: 142.8652
	trhc_S: -20.7711
	brhc_R: -115.8014
	brhc_A: -69.9204
	brhc_S: -57.0501
	  menc: 0
  normal_L: -106.3453
  normal_P: -134.3514
  normal_S: -26.3198
   ...
--> Data Axes Tilt:  Oblique (9.897 deg. from plumb)
--> Data Axes Approximate Orientation:
-->   first  (x) = Right-to-Left
-->   second (y) = Posterior-to-Anterior
-->   third  (z) = Inferior-to-Superior   [-orient RPI]
--> R-to-L extent:   -98.572 [R] -to-   114.428 [L] -step-     3.000 mm [ 72 voxels]
--> A-to-P extent:  -137.521 [A] -to-    75.479 [P] -step-     3.000 mm [ 72 voxels]
--> I-to-S extent:   -62.095 [I] -to-    87.902 [S] -step-     3.000 mm [ 51 voxels]
--> Geometry String: "MATRIX(2.99747,-0.109477,-0.056449,-98.57211,-0.098407,-2.955355,0.506141,75.47923,0.074081,0.503874,2.956386,-62.09487):72,72,51"
*/  
  save_nii(nii,strcat(pfile,'.nii'));

% If a {pfile}_tensor.dat file exists, generate {pfile}.bvec {PFILE}.bval files
  tensorfile = strcat(pfile,'_tensor.dat');
  if(exist(tensorfile, 'file') == 2)
    bvec = importdata(tensorfile);
    ndiff = bvec(1);
    bvec = reshape(bvec(2:end),[3,ndiff]);
    bvec(1,:) = -bvec(1,:);
    nT2 = size(d,4) - ndiff;
    bval = 1000*[zeros(1,nT2), ones(1,ndiff)]; % B1000 assumed.. find pfile flag (user1?)
    bvec = [zeros(3,nT2), bvec];
    fid = fopen([pfile,'.bval'],'w');
    fprintf(fid,'%.0f ',bval);
    fprintf(fid,'\n');
    fclose(fid);
    fid = fopen([pfile,'.bvec'],'w');
    for i=1:3,
      fprintf(fid,'%.3f ',bvec(i,:));
      fprintf(fid,'\n');
    end
    fclose(fid);
  end

  return;
  
/*  -------- incorrect matrix offsets ------------------------------------------>
$ 3dinfo P27648.7.nii.gz 
++ 3dinfo: AFNI version=AFNI_17.2.05 (Jul 27 2017) [64-bit]
** AFNI converts NIFTI_datatype=64 (FLOAT64) in file P27648.7.nii.gz to FLOAT32
     Warnings of this type will be muted for this session.
     Set AFNI_NIFTI_TYPE_WARN to YES to see them all, NO to see none.

Dataset File:    P27648.7.nii.gz
Identifier Code: NII_mcdQ3zsn2mMllBZDBDcUcw  Creation Date: Fri Jan 19 16:33:22 2018
Template Space:  ORIG
Dataset Type:    Echo Planar (-epan)
Byte Order:      LSB_FIRST {assumed} [this CPU native = LSB_FIRST]
Storage Mode:    NIFTI
Storage Space:   6,345,216 (6.3 million [mega]) bytes
Geometry String: "MATRIX(-3,0,0,-3,0,-3,0,-3,0,0,3,3):72,72,51"
Data Axes Tilt:  Plumb
Data Axes Orientation:
  first  (x) = Left-to-Right
  second (y) = Posterior-to-Anterior
  third  (z) = Inferior-to-Superior   [-orient LPI]
R-to-L extent:  -216.000 [R] -to-    -3.000 [R] -step-     3.000 mm [ 72 voxels]
A-to-P extent:  -216.000 [A] -to-    -3.000 [A] -step-     3.000 mm [ 72 voxels]
I-to-S extent:     3.000 [S] -to-   153.000 [S] -step-     3.000 mm [ 51 voxels]
Number of time steps = 6  Time step = 2.80000s  Origin = 0.00000s
  -- At sub-brick #0 '?' datum type is float
  -- At sub-brick #1 '?' datum type is float
  -- At sub-brick #2 '?' datum type is float
** For info on all 6 sub-bricks, use '3dinfo -verb' **

*/

/* -------------- With Matrix Offsets ---------------------------------------------> 

$ 3dinfo P27648.7_EpiMultiPhaseRecon.nii.gz 
++ 3dinfo: AFNI version=AFNI_17.2.05 (Jul 27 2017) [64-bit]

Dataset File:    P27648.7_EpiMultiPhaseRecon.nii.gz
Identifier Code: NII_7I9CsJ_xS0jkOMX3FA5XMw  Creation Date: Fri Jan 19 16:42:05 2018
Template Space:  ORIG
Dataset Type:    Echo Planar (-epan)
Byte Order:      LSB_FIRST {assumed} [this CPU native = LSB_FIRST]
Storage Mode:    NIFTI
Storage Space:   4,230,144 (4.2 million [mega]) bytes
Geometry String: "MATRIX(2.99747,-0.109477,-0.056449,-98.57211,-0.098407,-2.955355,0.506141,75.47923,0.074081,0.503874,2.956386,-62.09487):72,72,51"
Data Axes Tilt:  Oblique (9.897 deg. from plumb)
Data Axes Approximate Orientation:
  first  (x) = Right-to-Left
  second (y) = Posterior-to-Anterior
  third  (z) = Inferior-to-Superior   [-orient RPI]
R-to-L extent:   -98.572 [R] -to-   114.428 [L] -step-     3.000 mm [ 72 voxels]
A-to-P extent:  -137.521 [A] -to-    75.479 [P] -step-     3.000 mm [ 72 voxels]
I-to-S extent:   -62.095 [I] -to-    87.902 [S] -step-     3.000 mm [ 51 voxels]
Number of time steps = 8  Time step = 2.80000s  Origin = 0.00000s
  -- At sub-brick #0 '?' datum type is short
  -- At sub-brick #1 '?' datum type is short
  -- At sub-brick #2 '?' datum type is short
** For info on all 8 sub-bricks, use '3dinfo -verb' **

 */
 
 /*

>>pfile='P27648.7';
>> pfileHandle = GERecon('Pfile.Load', pfile);
>> header = GERecon('Pfile.Header', pfileHandle);
>> header

header = 

        RawHeader: [1x1 struct]
     DataAcqTable: [1x1 struct]
         ExamData: [1x1 struct]
       SeriesData: [1x1 struct]
        ImageData: [1x1 struct]
    PrescanHeader: [1x1 struct]
         GradData: [1x1 struct]
         CttEntry: [1x4 struct]

>> header.ImageData
 header.ImageData

ans = 

              autoSubParam: [1x1 struct]
            double_padding: [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]
                      dfov: 216
                 dfov_rect: 216
                    sctime: 28000008
                   slthick: 3
               scanspacing: 0
                       loc: -41.5775
                   tbldlta: 0
                       nex: 1
                   reptime: 0
                    saravg: 0.1356
                   sarpeak: 1.6031
                 pausetime: 0
                       vbw: 250
                     user0: 1
                     user1: 0
                     user2: 0
                     user3: 2
                     user4: 0
                     user5: 1
                     user6: 8004
                     user7: 1
                     user8: 1
                     user9: 0
                    user10: 0
                    user11: 1
                    user12: 0
                    user13: 0
                    user14: 0
                    user15: 0
                    user16: 0
                    user17: 2
                    user18: 4000
                    user19: 1
                    user20: 2
                    user21: 1
                    user22: 1
                  proj_ang: 0
                concat_sat: 0
                    user23: 0
                    user24: 0
                x_axis_rot: 0
                y_axis_rot: 0
                z_axis_rot: 0
                   ihtagfa: 0
                   ihtagor: 0
                   ihbspti: 0
                rtia_timer: 0
                       fps: 0
                 vencscale: 0
                      dbdt: 0
                   dbdtper: 100
                estdbdtper: 99.4493
                 estdbdtts: 77.8000
                saravghead: 0.6413
           neg_scanspacing: 0
                    user25: 0
                    user26: 0
                    user27: 0
                    user28: 0
                    user29: 0
                    user30: 0
                    user31: 0
                    user32: 0
                    user33: 0
                    user34: 0
                    user35: 0
                    user36: 0
                    user37: 0
                    user38: 0
                    user39: 0
                    user40: 0
                    user41: 0
                    user42: 0
                    user43: 0
                    user44: 0
                    user45: 0
                    user46: 0
                    user47: 0
                    user48: 0
              RegressorVal: 0
                SliceAsset: 0
                PhaseAsset: 0
                 sarValues: [0.1356 1.6031 0.6413 0]
                  shim_fov: [0 0]
                shim_ctr_R: [0 0]
                shim_ctr_A: [0 0]
                shim_ctr_S: [0 0]
                     dim_X: 72
                     dim_Y: 72
                 pixsize_X: 3
                 pixsize_Y: 3
                     ctr_R: -3.9513
                     ctr_A: 32.9297
                     ctr_S: -41.5775
                    norm_R: 0.0188
                    norm_A: -0.1687
                    norm_S: 0.9855
                    tlhc_R: 107.8988
                    tlhc_A: 135.7799
                    tlhc_S: -26.1049
                    trhc_R: -107.9191
                    trhc_A: 142.8652
                    trhc_S: -20.7711
                    brhc_R: -115.8014
                    brhc_A: -69.9204
                    brhc_S: -57.0501
                      menc: 0
                  normal_L: -106.3453
                  normal_P: -134.3514
                  normal_S: -26.3198
                       osf: 1
              fermi_radius: 36
               fermi_width: 10
                 fermi_ecc: 1
                   CSAccel: 0
                       esp: 0
                  bwPerPix: 0
                pixelSizeX: 0
                pixelSizeY: 0
            fov_freq_scale: 1
           fov_phase_scale: 1
             slthick_scale: 1
            freq_loc_shift: 0
           phase_loc_shift: 0
           slice_loc_shift: 0.4044
              b1PeakLowSAR: -9.9900
               b1rmsLowSAR: -9.9900
                  wbLowSAR: -9.9900
                headLowSAR: -9.9900
       seriesMaxTimeLowSAR: -9.9900
             float_padding: [0 0 0 0 0 0 0 0 0]
                cal_fldstr: 0
            user_usage_tag: 0
          user_fill_mapMSW: 0
          user_fill_mapLSW: 0
               im_archived: 0
               im_complete: 0
              int_padding1: [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]
               im_datetime: 1.5162e+09
              im_actual_dt: 1.5162e+09
                        tr: 2800000
                        ti: 0
                        te: 30000
                       te2: 0
                      tdel: 2854
                    mindat: 0
                   obplane: 18
                   slocfov: 0
                 obsolete1: 0
                 obsolete2: 0
               user_bitmap: 63308541
                      iopt: 67141632
              psd_datetime: 1.5156e+09
                 rawrunnum: 27648
                  intr_del: 0
                im_lastmod: 1.5162e+09
                  im_pds_a: 0
                  im_pds_c: 0
                  im_pds_u: 0
               thresh_min1: 0
               thresh_max1: 0
               thresh_min2: 0
               thresh_max2: 0
                  numslabs: 0
               locsperslab: 0
                  overlaps: 0
                slop_int_4: 0
                      dfax: 0
                    fphase: 8
                offsetfreq: 0
                   b_value: 0
                     iopt2: 0
                 ihtagging: 0
                  ihtagspc: 0
                 ihfcineim: 0
                 ihfcinent: 0
                   num_seg: 0
                   oprtarr: 0
                  averages: 0
             station_index: 0
             station_total: 0
                     iopt3: 0
                    delAcq: 0
                rxmbloblen: 1616
               rxmblob_pad: 0
                     im_no: 1
                     imgrx: 15
               temp_phases: 4
               driver_freq: 60
                driver_amp: 50
            driverCyc_Trig: 0
                   MEG_dir: 4
               rescan_time: 0
              spokesPerSeg: 0
              recoveryTime: 0
                  t2PrepTE: 0
                     hoecc: 0
                    rtb0cc: 0
              user_bitmap2: 0
            MultiBandAccel: 0
                  tissueT1: 0
            syn_acq_b_info: 0
              Acquirephase: 0
                  HKTAccel: 0
         saveOriginalImage: 0
      cardt1map_hb_pattern: 0
      distortionCorrection: 0
                ihdistcorr: 0
                 ihpepolar: 0
              int_padding2: [0 0 0 0 0 0 0 0 0]
                 imatrix_X: 72
                 imatrix_Y: 72
                   im_exno: 8421
                img_window: 0
                 img_level: 0
                   numecho: 1
                   echonum: 1
                   im_uniq: 0
                   im_seno: 6
                  contmode: 0
                     serrx: 1
              screenformat: 16
                     plane: 16
               im_compress: 1
              im_scouttype: 0
                    contig: 1
                   hrtrate: 0
                 trgwindow: 0
                   imgpcyc: 0
                 obsolete3: 0
                 obsolete4: 0
                 obsolete5: 0
                   mr_flip: 90
                    cphase: 0
                    swappf: 1
                  pauseint: 0
                 obsolete6: 0
                 obsolete7: 0
                 obsolete8: 0
                not_used_1: 0
                     imode: 1
                      pseq: 2
                  pseqmode: 0
                      ctyp: 3
                  surfctyp: 0
                  surfcext: 0
                 supp_tech: 0
                   slquant: 51
                      gpre: 0
                   satbits: 0
                      scic: 0
                  satxloc1: 9990
                  satxloc2: 9990
                  satyloc1: 9990
                  satyloc2: 9990
                  satzloc1: 9990
                  satzloc2: 9990
                 satxthick: 0
                 satythick: 0
                 satzthick: 0
                      flax: 0
                      venc: 0
               thk_disclmr: 0
                 obsolete9: 0
                obsolete10: 0
                image_type: 0
              vas_collapse: 0
                  proj_alg: 0
              echo_trn_len: 2
                 frac_echo: 2
                prep_pulse: 0
                 cphasenum: 0
                  var_echo: 0
                 scanactno: 1
                  vasflags: 2
                 integrity: 0
                  freq_dir: 1
                  vas_mode: 0
                   pscopts: 2
                obsolete11: 0
                obsolete12: 0
                obsolete13: 0
                unoriginal: 0
               interleaves: 2
              effechospace: 232
               viewsperseg: 0
                      rbpm: 0
                   rtpoint: 0
                  rcvrtype: 1
                   sarMode: 1
                  dBdtMode: 1
                   govBody: 1
             sarDefinition: 7
                no_shimvol: 0
             shim_vol_type: 0
             current_phase: 1
                 art_level: 0
        slice_group_number: 1
    number_of_slice_groups: 1
          show_in_autoview: 1
      slice_number_inGroup: 1
                   specnuc: 1
            label_duration: 0
            ihbsoffsetfreq: 4000
              scale_factor: 0
               volume_prop: 0
           excitation_mode: 0
             fat_shift_dir: 0
                ab1_status: 0
                    ab1_tg: -1
                      vest: 0
                  msde_dir: 0
           put_in_database: 1
                scout_scan: 0
             short_padding: [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]
                   psdname: 'muxarcepi                        '
                 proj_name: '             '
                 psd_iname: 'EPI          '
                 im_diskid: 0
                      pdid: '              '
                   im_suid: '3T  '
                contrastIV: '                 '
              contrastOral: '                 '
                   loc_ras: 73
                 forimgrev: '    '
                     cname: '32Ch Head        '
                im_verscre: '27'
                im_verscur: '27'
              im_alloc_key: 'MRIFCCALLOC  '
                   ref_img: 0
                   sum_img: 0
               filter_mode: '                '
                slop_str_2: '                '
                 image_uid: '                                '
                   sop_uid: '                                '
                   GEcname: 'C-GE_32Ch Head          '
              usedCoilData: '                                                                                …'
           astcalseriesuid: '+;"G*³µ¶94¶U²jg²brA;B    '
          purecalseriesuid: '+;"G*³µ¶94¶U²jg²brA;B    '
           xml_psc_shm_vol: '                                '
                   rxmpath: '                                                                '
              psdnameannot: '                                 '
                   anatomy: 'SRT%5CNone%5CT-D1100                                   '
           img_hdr_padding: '                                                                                …'

 
 */
 
