fprintf('\nMaking dprecor.mex and cprecor executables for the computation of precorrection matrices\n');
if ispc
    mex dprecor.cpp -largeArrayDims -v COMPFLAGS="/openmp /O2 /Wall $COMPFLAGS"
    mex cprecor.cpp -largeArrayDims -v COMPFLAGS="/openmp /O2 /Wall $COMPFLAGS"
end
if ismac
    mex dprecor.cpp -largeArrayDims -v CFLAGS="-fopenmp -Wall -O3 \$CFLAGS" LDFLAGS="-fopenmp \$LDFLAGS"
    mex cprecor.cpp -largeArrayDims -v CFLAGS="-fopenmp -Wall -O3 \$CFLAGS" LDFLAGS="-fopenmp \$LDFLAGS"
end
    