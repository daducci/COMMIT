classdef DICTIONARY_MakeOperator

properties
	nS					% number of SAMPLES
	nF					% number of FIBERS
	nR					% number of FIBER RADII
	nE					% number of EC compartments
	nT					% number of EC TORTUOSITY values
	nV					% number of VOXELS
	nI					% number of ISO compartments
    adjoint
end


methods

	function obj = DICTIONARY_MakeOperator()
		global DICTIONARY KERNELS
		obj.adjoint = 0;
        obj.nS = KERNELS.nS;
        obj.nF = DICTIONARY.IC.nF;
		obj.nR = numel(KERNELS.wmr);
		obj.nE = DICTIONARY.EC.nE;
		obj.nT = numel(KERNELS.wmh);
		obj.nV = DICTIONARY.nV;
		obj.nI = numel(KERNELS.iso);
    end


    function obj = ctranspose( obj )
        obj.adjoint = 1-obj.adjoint;
	end


    function res = size( obj, dim )
		res = [obj.nV*obj.nS obj.nR*obj.nF+obj.nT*obj.nE+obj.nI*obj.nV];
		if obj.adjoint == 1
			res = res([2 1]);
		end
		if nargin==2
			res = res( dim );
		end
	end


    function v = mtimes( obj, x )
		global DICTIONARY KERNELS THREADS
		if ~isa(x,'double'), x = double(x); end
        if obj.adjoint == 0
			% A*x
			if ( obj.size(2) ~= size(x,1) ), error('A.times(): dimensions do not match'); end
			v = DICTIONARY_MakeOperator_A( DICTIONARY, KERNELS, x, THREADS );
		else
			% At*y
			if ( obj.size(2) ~= size(x,1) ), error('At.times(): dimensions do not match'); end
			v = DICTIONARY_MakeOperator_At( DICTIONARY, KERNELS, x, THREADS );
        end
	end

end

end
