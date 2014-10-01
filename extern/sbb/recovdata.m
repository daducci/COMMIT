function data = recovdata(out)
% function data = recovdata(out)

    ds = RandStream.getDefaultStream;
    ds.State = out.seed{1};
    if out.sp1 == 1
        data.A = rand(out.m, out.n);
    else
        data.A = sprand(out.m, out.n, out.sp1);
    end

    ds.State = out.seed{2};
    data.xt = sprand(out.n, 1, out.sp2);
    b1 = data.A * data.xt;
    
    idx = find(data.xt == 0);
    y = zeros(size(data.xt));

    ds.State = out.seed{3};
    y(idx) = .01 * rand(length(idx), 1);

    t = data.A' * b1 - y;
    
    if issparse(data.A)
        R = qr(data.A');
    else
        R = triu(qr(data.A'));
    end
    data.b = R\(R'\(data.A*t));
    r = t - data.A'*data.b;
    e = R\(R'\(data.A*r));
    data.b = data.b + e;

