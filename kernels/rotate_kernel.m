function [ Kr ] = rotate_kernel( filename, scheme, AUX, idx_IN, idx_OUT, Ylm_OUT )
	fid = fopen( filename, 'r', 'b' );
	K = fread(fid,'float');
	fclose(fid);

	% rotate it
	Kr = zeros( size(scheme.camino,1), 181, 181, 'single' );
	for ox = 1:181
	for oy = 1:181
		Ylm_rot = AUX.Ylm_rot{ ox, oy };
		for s = 1:numel(scheme.shells)
			Klm = AUX.fit * K( idx_IN{s} );		% fti SH to shell to rotate
			Rlm = zeros( size(Klm) );
			idx = 1;
			for l = 0 : 2 : AUX.lmax
				const = sqrt(4.0*pi/(2.0*l+1.0)) * Klm( (l*l + l + 2.0)/2.0 );
				for m = -l : l
					Rlm(idx) = const * Ylm_rot(idx);
					idx = idx+1;
				end
			end
			Kr( idx_OUT{s}, ox, oy ) = Ylm_OUT{s} * Rlm;
		end
	end
	end
