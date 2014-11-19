fprintf('\n-> Processing kernels...\n');
ticID = tic;

% Ensure double data type for mex code
% ------------------------------------
fprintf('\t - convert to double data type...');
for i = 1:numel(KERNELS.wmr)
	KERNELS.wmr{i} = double( KERNELS.wmr{i}(CONFIG.scheme.dwi_idx,:,:) );
end
for i = 1:numel(KERNELS.wmh)
	KERNELS.wmh{i} = double( KERNELS.wmh{i}(CONFIG.scheme.dwi_idx,:,:) );
end
for i = 1:numel(KERNELS.iso)
	KERNELS.iso{i} = double( KERNELS.iso{i}(CONFIG.scheme.dwi_idx,:,:) );
end
KERNELS.nS = CONFIG.scheme.dwi_count;
fprintf( ' [ OK ]\n' );


% De-mean kernels as LiFE
% -----------------------
if ( CONFIG.doDemean )
	fprintf('\t - de-meaning atoms as LiFE...');
	for j = 1:181
	for k = 1:181
		for i = 1:numel(KERNELS.wmr)
			KERNELS.wmr{i}(:,j,k) = KERNELS.wmr{i}(:,j,k) - mean(KERNELS.wmr{i}(:,j,k),1);
		end
		for i = 1:numel(KERNELS.wmh)
			KERNELS.wmh{i}(:,j,k) = KERNELS.wmh{i}(:,j,k) - mean(KERNELS.wmh{i}(:,j,k),1);
		end
		for i = 1:numel(KERNELS.iso)
			KERNELS.iso{i}(:,j,k) = KERNELS.iso{i}(:,j,k) - mean(KERNELS.iso{i}(:,j,k),1);
		end
	end
	end
	fprintf( ' [ OK ]\n' );
end


% Add b0 reference (promote unitary sum)
% --------------------------------------
if ( CONFIG.useReference )
	fprintf('\t - adding reference...');
	for i = 1:numel(KERNELS.wmr)
		KERNELS.wmr{i} = KERNELS.wmr{i}([1 1:KERNELS.nS],:,:);
		KERNELS.wmr{i}(1,:,:) = 1;
	end
	for i = 1:numel(KERNELS.wmh)
		KERNELS.wmh{i} = KERNELS.wmh{i}([1 1:KERNELS.nS],:,:);
		KERNELS.wmh{i}(1,:,:) = 1;
	end
	for i = 1:numel(KERNELS.iso)
		KERNELS.iso{i} = KERNELS.iso{i}([1 1:KERNELS.nS],:,:);
		KERNELS.iso{i}(1,:,:) = 1;
	end
	KERNELS.nS = KERNELS.nS + 1;
	fprintf( ' [ OK ]\n' );
end


% Normalize atoms
% ---------------
if ( CONFIG.kernels.doNormalize )
	fprintf('\t - normalizing to have unitary norm...');
	for i = 1:numel(KERNELS.wmr)
		KERNELS.wmr_norm(i) = norm( KERNELS.wmr{i}(:,1,1) );
		for j = 1:181
		for k = 1:181
			KERNELS.wmr{i}(:,j,k) = KERNELS.wmr{i}(:,j,k) / KERNELS.wmr_norm(i);
		end
		end
	end
	for i = 1:numel(KERNELS.wmh)
		KERNELS.wmh_norm(i) = norm( KERNELS.wmh{i}(:,1,1) );
		for j = 1:181
		for k = 1:181
			KERNELS.wmh{i}(:,j,k) = KERNELS.wmh{i}(:,j,k) / KERNELS.wmh_norm(i);
		end
		end
	end
	for i = 1:numel(KERNELS.iso)
		KERNELS.iso_norm(i) = norm( KERNELS.iso{i}(:,1,1) );
		KERNELS.iso{i} = KERNELS.iso{i} / KERNELS.iso_norm(i);
	end
	fprintf( ' [ OK ]\n' );
end


fprintf( '   [ %.2f seconds ]\n', toc(ticID) );
clear ticID i j k
