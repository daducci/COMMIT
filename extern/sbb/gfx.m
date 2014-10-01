function [g flag] = gfx(x)
    global data

    flag = 0;
    Ax = data.A * x;
    g = (data.A' * (Ax - data.b));

