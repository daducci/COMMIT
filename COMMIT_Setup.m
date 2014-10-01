%
% Initialization for COMMIT
%

% Global variables
% ================
global COMMIT_path CAMINO_path
global KERNELS DICTIONARY THREADS


% Path definition
% ===============
COMMIT_path = '/Users/Ale/Documents/LTS5/Projects/Tractography/COMMIT';
CAMINO_path = '/Users/Ale/Downloads/camino/bin/';

addpath( fullfile(COMMIT_path,'code','kernels') )
addpath( fullfile(COMMIT_path,'code','dictionary') )
addpath( fullfile(COMMIT_path,'code','optimization') )
addpath( fullfile(COMMIT_path,'code','extern','sbb') )
addpath( fullfile(COMMIT_path,'code','extern','spgl1') )
