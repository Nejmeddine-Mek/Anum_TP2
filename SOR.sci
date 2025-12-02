    // fix this tthing samsd
    
function [x,ops, iter] = SOR(A,b,epsilon)
    // first of all, we need to verify that A is symmetric and positive definite
    isSymmetric = norm(A - A',"fro") < 1e-12
    
    isPosDef = %t
    try
        chol(A)
    catch
        isPosDef = %f
    end
    if isSymmetric & isPosDef then
        disp("can apply sor")

    else
        disp("can''t apply sor")
        return
    end
    
    D = diag(diag(A))
    L = tril(A,-1)
    U = triu(A,1)
    
    // at this point, one shall start executing the SOR method as follows
    // initially set, the parameter w to 1
    w = 1
    // and then we will use the residual and 
    R = []
    W = [1]
    X_cur = zeros(length(b),1)
    r = b - A * X_cur
    R($+1) = norm(r)
    n = size(A,1);
    //Lsor = inv(1 ./ w .* D - L) * (U - (1 - 1 ./ w) * D)
   // X_next = Lsor * X_cur + inv(1 ./ w .* D - L) * b
    while norm(r) >= epsilon do
           
        X_old = X_cur;
        
                // --- SOR update ---
        for i = 1:n
            // left sum
            if i > 1 then
                sigma1 = A(i,1:i-1) * X_cur(1:i-1);
            else
                sigma1 = 0;
            end
        
            // right sum
            if i < n then
                sigma2 = A(i,i+1:n) * X_old(i+1:n);
            else
                sigma2 = 0;
            end
        
            X_cur(i) = (1-w)*X_old(i) + w * (b(i) - sigma1 - sigma2) / A(i,i);
        end

        
        r_new = b - A * X_cur;
        
        R($+1) = norm(r_new);
        
        if norm(r_new) < norm(r) then
            w = min(w + 0.1, 1.9);
            W($+1) = w;
        end
        
        r = r_new;
    end

   // now plotting R
    figure()
    plot(R)
    xtitle("Residual Norm vs Iteration", "Iteration", "||r||")
    figure()
    plot(W)
    xtitle("parametr w vs Iteration", "Iteration", "||r||")
    
    iter = W
    ops = R
    x = X_cur
endfunction
[x,ops,iter] = SOR([4 1 1; 1 3 -1; 1 -1 2],[6;3;2],0.001)

disp(x, ops, iter)
