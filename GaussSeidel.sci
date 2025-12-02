function [x,ops,iter] = GaussSeidel(A,b,epsilon)
    D = diag(diag(A))
    L = -1 .* tril(A,-1)
    U = -1 .* triu(A,1)
    Lgs = inv(D - L)
    if max(abs(spec(Lgs))) >= 1 then
        x = 0
        ops = 0
        iter = 0
        return
    end
    X_cur = zeros(length(b),1)
    X_next = Lgs * U * X_cur + Lgs * b
    while norm(X_next - X_cur) >= epsilon do
        X_cur = X_next
        X_next = Lgs * U * X_cur + Lgs * b
    end
    ops = 1
    iter = 1
    x = X_next
endfunction

    [x,y,z] = GaussSeidel([4 1 1; 1 3 -1; 1 -1 2],[6;3;2],0.001)
    disp(x)
