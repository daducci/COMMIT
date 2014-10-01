function gensyn()
% GENSYN()
% Generate synthetic NNLS datasets
%
% rs.dense : info. on dense problems
% rs.sparse : info. on sparse problems
%
% rs.*.act : number of active variables at the solution
% rs.*.val : inf-norm of the projected gradient at the solution
% rs.*.out : random number info to reproduce datasets
%

    datadir = '../../data/';

    rs.dense.act = zeros(6,1);
    rs.dense.val = zeros(6,1);
    rs.dense.out = cell(6,1);

    drow = 300;
    dcol = 200;

    for i = 1 : 6
        fname = sprintf([datadir 'syn_%d.mat'], i);
        [data rs.dense.out{i}] = gennnls(drow * 2^i, dcol * 2^i, 1, 0.3);
        rs.dense.act(i) = length(data.xt) - nnz(data.xt);
        rs.dense.val(i) = pgtol(data.A, data.b, data.xt);
        save(fname, 'data');
        clear data;
        fprintf('*');
    end

    rs.sparse.act = zeros(6,1);
    rs.sparse.val = zeros(6,1);
    rs.sparse.out = cell(6,1);

    srow = 25600;
    scol = 9600;

    for i = 1 : 6
        fname = sprintf([datadir 'syn_%d.mat'], 10 + i);
        [data rs.sparse.out{i}] = gennnls(srow, scol, 0.005*i, 0.3);
        rs.sparse.act(i) = length(data.xt) - nnz(data.xt);
        rs.sparse.val(i) = pgtol(data.A, data.b, data.xt);
        save(fname, 'data');
        clear data;
        fprintf('-');
    end
        
    save([datadir 'syn_rs.mat'], 'rs');

function [data out] = gennnls(m, n, sp1, sp2)
% function [data out] = gennnls(m, n, sp1, sp2)
% m : number of rows in A,
% n : number of cols in A,
% sp1 : sparsity of A,
% sp2 : sparsity of xt (true x),
%
    out.m = m;
    out.n = n;
    out.sp1 = sp1;
    out.sp2 = sp2;
    out.seed = cell(3,1);

    ds = RandStream.getDefaultStream;
    out.seed{1} = ds.State;
    if sp1 == 1
        data.A = rand(m, n);
    else
        data.A = sprand(m, n, sp1);
    end

    out.seed{2} = ds.State;
    data.xt = sprand(n, 1, sp2);
    data.bd = data.A * data.xt;
    
    idx = find(data.xt == 0);
    y = zeros(size(data.xt));

    out.seed{3} = ds.State;
    y(idx) = .01 * rand(length(idx), 1);

    t = data.A' * data.bd - y;
    
    if issparse(data.A)
        R = qr(data.A');
    else
        R = triu(qr(data.A'));
    end
    data.b = R\(R'\(data.A*t));
    r = t - data.A'*data.b;
    e = R\(R'\(data.A*r));
    data.b = data.b + e;


