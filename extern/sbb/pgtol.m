function val = pgtol(A, b, x)

grad = A' * (A * x - b);
idx = find(x ~= 0 | (x == 0 & grad < 0));
val = norm(abs(grad(idx)), 'inf');
    
