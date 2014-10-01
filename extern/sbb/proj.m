function [x flag] = proj(x)
flag = 0;
x(x<0) = 0;
%x(x < 1e-7) = 0;
    
