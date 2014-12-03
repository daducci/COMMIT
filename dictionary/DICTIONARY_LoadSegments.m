ticID = tic;
fprintf( '\n-> Loading the dictionary\n   ======================\n' );

DICTIONARY = [];
DICTIONARY.trkFilename = 'dictionary_fibers.trk';
DICTIONARY.dim    = niiSIGNAL.hdr.dime.dim(3:5);
DICTIONARY.pixdim = niiSIGNAL.hdr.dime.pixdim(3:5);

niiMASK = load_untouch_nii( fullfile(CONFIG.TRACKING_path,'dictionary_mask.nii') );
DICTIONARY.MASK = niiMASK.img;
clear niiMASK

% intra-axonal compartments
% -------------------------
fprintf( '\t- intra-axonal compartments...' );

DICTIONARY.IC = [];

pDict_IC_trkLen = fopen( fullfile(CONFIG.TRACKING_path,'dictionary_IC_trkLen.dict'), 'rb' );
if pDict_IC_trkLen < 0
	error( 'Dictionary not found. Launch ''trk2dictionary'' script first!' );
end
pDict_IC_f      = fopen( fullfile(CONFIG.TRACKING_path,'dictionary_IC_f.dict'),      'rb' );
pDict_IC_vx     = fopen( fullfile(CONFIG.TRACKING_path,'dictionary_IC_vx.dict'),     'rb' );
pDict_IC_vy     = fopen( fullfile(CONFIG.TRACKING_path,'dictionary_IC_vy.dict'),     'rb' );
pDict_IC_vz     = fopen( fullfile(CONFIG.TRACKING_path,'dictionary_IC_vz.dict'),     'rb' );
pDict_IC_ox     = fopen( fullfile(CONFIG.TRACKING_path,'dictionary_IC_ox.dict'),     'rb' );
pDict_IC_oy     = fopen( fullfile(CONFIG.TRACKING_path,'dictionary_IC_oy.dict'),     'rb' );
pDict_IC_len    = fopen( fullfile(CONFIG.TRACKING_path,'dictionary_IC_len.dict'),    'rb' );

fseek(pDict_IC_trkLen,0,'eof');
DICTIONARY.IC.nF = ftell(pDict_IC_trkLen) / 4;
fseek(pDict_IC_trkLen,0,'bof');

fseek(pDict_IC_f,0,'eof');
DICTIONARY.IC.n = ftell(pDict_IC_f) / 4;
fseek(pDict_IC_f,0,'bof');

DICTIONARY.IC.trkLen = fread( pDict_IC_trkLen, [DICTIONARY.IC.nF 1], '*float32'   );
DICTIONARY.IC.fiber  = fread( pDict_IC_f,      [DICTIONARY.IC.n  1], '*uint32'  );	% fiber containing this segments

vx     = fread( pDict_IC_vx,     [DICTIONARY.IC.n  1], '*uint8'   );	% position of the segment
vy     = fread( pDict_IC_vy,     [DICTIONARY.IC.n  1], '*uint8'   );
vz     = fread( pDict_IC_vz,     [DICTIONARY.IC.n  1], '*uint8'   );
DICTIONARY.IC.v  = uint32(vx)  + DICTIONARY.dim(1) * ( uint32(vy)  + DICTIONARY.dim(2) * uint32(vz) );
clear vx vy vz

ox     = fread( pDict_IC_ox,     [DICTIONARY.IC.n  1], '*uint8'   );	% orientation of the segment
oy     = fread( pDict_IC_oy,     [DICTIONARY.IC.n  1], '*uint8'   );
DICTIONARY.IC.o  = uint16(ox) + 181*uint16(oy);
clear ox oy

DICTIONARY.IC.len    = fread( pDict_IC_len,    [DICTIONARY.IC.n  1], '*float32' );	% length of the segment

if ( CONFIG.kernels.doNormalize )
	% all the columns will have same length
	for s = 1:DICTIONARY.IC.n
		f = DICTIONARY.IC.fiber(s) + 1;
		DICTIONARY.IC.len(s) = DICTIONARY.IC.len(s) / DICTIONARY.IC.trkLen(f);
	end
end

fclose(pDict_IC_trkLen);
fclose(pDict_IC_f);
fclose(pDict_IC_vx);
fclose(pDict_IC_vy);
fclose(pDict_IC_vz);
fclose(pDict_IC_ox);
fclose(pDict_IC_oy);
fclose(pDict_IC_len);
clear pDict_IC_*

% reorder the segments based on the "v" field
[DICTIONARY.IC.v, idx] = sort( DICTIONARY.IC.v );
DICTIONARY.IC.o     = DICTIONARY.IC.o( idx );
DICTIONARY.IC.fiber = DICTIONARY.IC.fiber( idx );
DICTIONARY.IC.len   = DICTIONARY.IC.len( idx );
clear idx

fprintf( '\t[ %d fibers and %d segments ]\n', DICTIONARY.IC.nF, DICTIONARY.IC.n );


% extra-axonal compartments
% -------------------------
fprintf( '\t- extra-axonal compartments...' );

DICTIONARY.EC = [];

pDict_EC_vx  = fopen( fullfile(CONFIG.TRACKING_path,'dictionary_EC_vx.dict'),     'rb' );
pDict_EC_vy  = fopen( fullfile(CONFIG.TRACKING_path,'dictionary_EC_vy.dict'),     'rb' );
pDict_EC_vz  = fopen( fullfile(CONFIG.TRACKING_path,'dictionary_EC_vz.dict'),     'rb' );
pDict_EC_ox  = fopen( fullfile(CONFIG.TRACKING_path,'dictionary_EC_ox.dict'),     'rb' );
pDict_EC_oy  = fopen( fullfile(CONFIG.TRACKING_path,'dictionary_EC_oy.dict'),     'rb' );

fseek(pDict_EC_vx,0,'eof');
DICTIONARY.EC.nE = ftell(pDict_EC_vx);
fseek(pDict_EC_vx,0,'bof');

vx     = fread( pDict_EC_vx,  [DICTIONARY.EC.nE  1], '*uint8'   );	% position of the segment
vy     = fread( pDict_EC_vy,  [DICTIONARY.EC.nE  1], '*uint8'   );
vz     = fread( pDict_EC_vz,  [DICTIONARY.EC.nE  1], '*uint8'   );
DICTIONARY.EC.v  = uint32(vx)  + DICTIONARY.dim(1) * ( uint32(vy)  + DICTIONARY.dim(2) * uint32(vz) );
clear vx vy vz

ox     = fread( pDict_EC_ox,     [DICTIONARY.EC.nE  1], '*uint8'   );	% orientation of the segment
oy     = fread( pDict_EC_oy,     [DICTIONARY.EC.nE  1], '*uint8'   );
DICTIONARY.EC.o  = uint16(ox) + 181*uint16(oy);
clear ox oy

fclose(pDict_EC_vx);
fclose(pDict_EC_vy);
fclose(pDict_EC_vz);
fclose(pDict_EC_ox);
fclose(pDict_EC_oy);
clear pDict_EC_*

% reorder the segments based on the "v" field
[DICTIONARY.EC.v, idx] = sort( DICTIONARY.EC.v );
DICTIONARY.EC.o     = DICTIONARY.EC.o( idx );
clear idx

fprintf( '\t[ %d segments ]\n', DICTIONARY.EC.nE );


% isotropic compartment
% ---------------------
fprintf( '\t- isotropic compartments...' );

DICTIONARY.nV = nnz( DICTIONARY.MASK );
[ vx, vy, vz ] = ind2sub( niiSIGNAL.hdr.dime.dim(3:5), find( DICTIONARY.MASK ) );
DICTIONARY.ISO.v = uint32(vx-1) + DICTIONARY.dim(1) * ( uint32(vy-1) + DICTIONARY.dim(2) * uint32(vz-1) );

% reorder the segments based on the "v" field
DICTIONARY.ISO.v = sort( DICTIONARY.ISO.v );

clear vx vy vz

fprintf( '\t[ %d voxels ]\n', DICTIONARY.nV );


% post-processing
% ---------------
fprintf( '\t- post-processing...' );

DICTIONARY.MASKidx = find( permute(repmat(DICTIONARY.MASK,[1 1 1 niiSIGNAL.hdr.dime.dim(2) ]),[4 1 2 3]) );

idx = find( DICTIONARY.MASK );
lut = zeros( DICTIONARY.dim, 'uint32' );
for i = 1:numel(idx)
    lut( idx(i) ) = (i-1)*KERNELS.nS;
end


for i = 1:numel(DICTIONARY.IC.v)
    DICTIONARY.IC.v(i) = lut( DICTIONARY.IC.v(i)+1 );
end
for i = 1:numel(DICTIONARY.EC.v)
    DICTIONARY.EC.v(i) = lut( DICTIONARY.EC.v(i)+1 );
end
for i = 1:numel(DICTIONARY.ISO.v)
    DICTIONARY.ISO.v(i) = lut( DICTIONARY.ISO.v(i)+1 );
end

clear idx lut i

fprintf( '\t\t[ OK ]\n' );


fprintf( '   [ %.2f seconds ]\n', toc(ticID) );
clear ticID
