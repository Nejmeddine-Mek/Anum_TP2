function [x,ops,iter] = Jacobi(A,b,epsilon)
    D = diag(diag(A))
    L = -1 .* tril(A,-1)
    U = -1 .* triu(A,1)
    D_inv = diag( 1 ./ diag(D))
    Lj = D_inv * ( L + U)
    radius = max(abs(spec(Lj)))
    if radius < 1 then
        disp("ok")
    else
        disp("no")
        // one quits here
    end
    
    X_cur = zeros(length(b),1)
    X_next = Lj * X_cur + D_inv * b
    ops = 3
    iter = 1
    while norm(X_next - X_cur) >= epsilon do
        X_cur = X_next
        X_next = Lj * X_cur + D_inv * b
        ops = ops + 3
        iter = iter + 1
    end
    
    x = X_next
endfunction

[x,ops,iter] = Jacobi([4 1 1; 1 3 -1; 1 -1 2],[6;3;2],0.001)
disp("result: ",x)
disp("iterations: ",iter)
disp("operations: ",ops)
