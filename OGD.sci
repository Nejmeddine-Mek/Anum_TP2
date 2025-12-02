
function x = OGD(A,b,epsilon)
    X = zeros(length(b),1);
    r = b - A * X;
    while(norm(r,2) >= epsilon) do
        a = (r' * r) / ((A *r)' * r);
        X = X + a * r;
        r = b - A * X;
    end
    x = X;
endfunction
