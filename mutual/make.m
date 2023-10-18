fprintf('\nMaking mutual.mex executable for the computation of near-by 1/r interactions\n');
if ispc
    mex main.c cuber.c mutual.c -v -output mutual COMPFLAGS="/openmp /O2 /Wall $COMPFLAGS"
end
if ismac || isunix
    mex main.c cuber.c mutual.c -v -output mutual CFLAGS="-fopenmp -Wall -O3 \$CFLAGS" LDFLAGS="-fopenmp \$LDFLAGS"
end
    