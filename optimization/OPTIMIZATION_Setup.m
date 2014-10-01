ticID = tic;
fprintf( '\n-> Creating the linear operator\n   ==========================\n' );


% distribute DIFFERENT VOXELS to DIFFERENT THREAD to avoid conflicts
% ==================================================================
THREADS = (CONFIG.OPTIMIZATION.nTHREADS-1)*ones( DICTIONARY.IC.n, 1, 'uint8' );
if CONFIG.OPTIMIZATION.nTHREADS>1
	
	fprintf( '\t- distributing voxels to %d different threads...\n', CONFIG.OPTIMIZATION.nTHREADS );
	[Vsort Vidx] = sort( DICTIONARY.IC.v );
	[~,~,Uidx] = unique(Vsort,'stable');
	clear Vsort
	C = accumarray(Uidx,1);
	tID = 0;
	tot = 0;
	iPrev = 1;
	for i = 1:numel(C)
		if tot >= floor(DICTIONARY.IC.n/CONFIG.OPTIMIZATION.nTHREADS)
			THREADS( Vidx(Uidx>=iPrev & Uidx<i) ) = tID;
			tID   = tID+1;
			if tID==CONFIG.OPTIMIZATION.nTHREADS-1, break, end
			iPrev = i;
			tot   = C(i);
		else
			tot = tot + C(i);
		end
	end
end
clear Vidx Uidx tID tot iPrev i C


% Compiling the mex-files encoding the linear operator A
% ======================================================
fprintf( '\t- compiling mex files...\n' );
clear A At

mex( ...
	sprintf('-DnIC=%d',numel(KERNELS.wmr)), sprintf('-DnEC=%d',numel(KERNELS.wmh)), sprintf('-DnISO=%d',numel(KERNELS.iso)), ...
	sprintf('-DnTHREADS=%d',CONFIG.OPTIMIZATION.nTHREADS), '-O', '-silent', ...
	'-outdir', fullfile(COMMIT_path,'code','dictionary'), fullfile(COMMIT_path,'code','dictionary','DICTIONARY_MakeOperator_A.cpp') ...
);

mex( ...
	sprintf('-DnIC=%d',numel(KERNELS.wmr)), sprintf('-DnEC=%d',numel(KERNELS.wmh)), sprintf('-DnISO=%d',numel(KERNELS.iso)), ...
	sprintf('-DnTHREADS=%d',CONFIG.OPTIMIZATION.nTHREADS), '-O', '-silent', ...
	'-outdir', fullfile(COMMIT_path,'code','dictionary'), fullfile(COMMIT_path,'code','dictionary','DICTIONARY_MakeOperator_At.cpp') ...
);

A  = DICTIONARY_MakeOperator();
At = A';

Y = niiSIGNAL.img( DICTIONARY.MASKidx );

fprintf( '   [ %.2f seconds ]\n', toc(ticID) );
clear ticID
