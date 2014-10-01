function [f flag] = fx(x)
    global data
    global startTime

    flag = 0;

    Ax = data.A * x;
    f = 0.5 * norm(Ax - data.b)^2;
