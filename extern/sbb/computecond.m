function computecond()
% COMPUTECOND()
% Compute condition numbers 
%
% cs.dense : info. on dense problems
% cs.sparse : info. on sparse problems
%

    datadir = '../../data/';

    cs.odense = zeros(6,1);
    cs.rdense = zeros(6,1);
    cs.osparse = zeros(6,1);
    cs.rsparse = zeros(6,1);

    for i = [1 : 6]
        fname = sprintf([datadir 'syn_%d.mat'], i);
        load(fname);
        H = data.A' * data.A;
        cs.odense(i) = cond(H);
        idx = find(data.xt ~= 0);
        cs.rdense(i) = cond(H(idx, idx));
        fprintf('ocond[%d] = %f, rcond[%d] = %f\n', i, cs.odense(i), ...
                i, cs.rdense(i));
    end

    for i = 1 : 6 
        fname = sprintf([datadir 'syn_%d.mat'], 10 + i);
        load(fname);
        H = data.A' * data.A;
        cs.osparse(i) = condest(H);
        idx = find(data.xt ~= 0);
        cs.rsparse(i) = condest(H(idx, idx));
        fprintf('ocond[%d] = %f, rcond[%d] = %f\n', 10 + i, cs.osparse(i), ...
                10 + i, cs.rsparse(i));
    end

    save('syn_cs.mat', 'cs');
        
