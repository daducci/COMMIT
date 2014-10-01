
fprintf('build ASA\n');
cd ../ASA_CG-2.2/
mex -largeArrayDims asamex.c asa_cg.c

fprintf('build lbfgs-b\n');
%cd ../lbfgsb-for-matlab/
%mex -largeArrayDims -lstdc++ -o lbfgsb arrayofmatrices.cpp lbfgsb.cpp matlabexception.cpp matlabmatrix.cpp matlabprogram.cpp matlabscalar.cpp matlabstring.cpp program.cpp solver.f

cd ../lbfgs-1.1
mex -largeArrayDims lbfgs_mex.c list.c utils.c routines.f

fprintf('build SPG\n');

cd ../spg
mex FC=g95 -largeArrayDims spgmmex.F spg.f


cd ../common


fprintf('Done.\n');
