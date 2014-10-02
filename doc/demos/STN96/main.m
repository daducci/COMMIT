%
% DEMO: comparison between COMMIT and LIFE on STN96 data from the LiFE original pubblication
%
% Download data for the demo:
%   1) Download precomputed data (CSD reconstruction + probabilistic tracking) from:
%      http://hardi.epfl.ch/static/data/COMMIT_demos/STN96_scan1.zip
%   2) Unzip the file into the "STN96/scan1" folder in your data directory, e.g.
%      data/STN96/scan1/
%   3) Download the original dwi data from:
%      https://stacks.stanford.edu/file/druid:cs392kv3054/life_demo_data.tar.gz
%   4) Unzip the file "life_demo_scan1_subject1_b2000_150dirs_stanford.nii.gz"
%      and move it into the "scan1" folder with the name DWI.nii, i.e. 
%      gunzip life_demo_scan1_subject1_b2000_150dirs_stanford.nii.gz
%      mv life_demo_scan1_subject1_b2000_150dirs_stanford.nii data/STN96/scan1/DWI.nii


%% STEP 1: process the data with LiFE
step1_LIFE


%% STEP 2: process the data with COMMIT
step2_COMMIT


%% STEP 3: compare the two models
step3_compare
