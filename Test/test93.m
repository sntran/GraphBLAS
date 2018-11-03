function test93
%TEST93 test dpagerank and ipagerank

rng ('default') ;
addpath ('../Demo/MATLAB') ;

for n = [10 100 1000 1e4 1e5]

    fprintf ('\n--------------n: %d\n', n) ;
    nz = 8*n ;
    d = nz / n^2 ;
    A = sprand (n, n, d) ;
    A = spones (A) ;

    tic
    [r1, i1] = dpagerank (A) ;
    t1 = toc ;

    [r2, i2] = GB_mex_dpagerank (A) ;
    t2 = gbresults ;

    % [i1' i2']

    % results won't be identical because of different random number generators
    mismatch = length (find (i1 ~= i2)) ;
    e = norm (r1 - r2) / norm (r1) ;
    fprintf ('i1 i2 mismatch: %d\n', mismatch) ;
    fprintf ('r1-r2 = %g\n', e) ;
    fprintf ('time: MATLAB %g GraphBLAS %g speedup %g\n', t1, t2, t1/t2) ;
    assert (mismatch < 100) ;
    assert (e < 1e-6) ;

    fprintf ('\ninteger version:\n') ;
    
    % test the integer versions
    tic
    [ir1, ii1] = ipagerank (A) ;
    ti1 = toc ;

    [ir2, ii2] = GB_mex_ipagerank (A) ;
    ti2 = gbresults ;

    % [ii1' ii2']

    ir1 = ir1 / norm (ir1) ;
    ir2 = ir2 / norm (ir2) ;

    % results won't be identical because of different random number generators
    mismatch = length (find (ii1 ~= ii2)) ;
    e = norm (ir1 - ir2) / norm (ir1) ;
    fprintf ('i1 i2 mismatch: %d\n', mismatch) ;
    fprintf ('r1-r2 = %g\n', e) ;
    fprintf ('time: MATLAB %g GraphBLAS %g speedup %g\n', ti1, ti2, ti1/ti2) ;
    if (n < 1e4)
        assert (mismatch < n/10) ;
        assert (e < 1e-4) ;
    end

end

fprintf ('test93: all tests passed\n') ;

