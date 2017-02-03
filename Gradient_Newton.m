function [] = newtonDescent()
    n = 100;
    m = 200;
    randn('state',1);
    A=randn(m,n);
    alpha = .01;
    beta = .5;
    f_vals = [];
    iter_vals = [];
    eta_vals = [];
    step_vals = [];
    eta = 0.00001;
    
    x = randn(n,1)*.001;
    while(~dom(x,A))
        x = randn(n,1)*.001;
    end
    grad_x = grad(x,A);
    hess_x_inv = pinv(hess(x,A));
    maxIter = 1000;
    iter = 1;
    lambdaSquare = grad_x'*hess_x_inv*grad_x;
    while(lambdaSquare/2 > eta && iter<maxIter)
        direction = -hess_x_inv*grad_x;
        t = 1;
        iter_vals = [iter_vals;iter];
        f_vals = [f_vals;f(x,A)-(-144.6979)];
        eta_vals = [eta_vals;lambdaSquare];
        
        while(~dom(x+t*direction,A))
                t = beta*t;
        end
        while(f(x+t*direction,A) > f(x,A) + alpha*t*grad_x'*direction)
            t = beta*t;
        end
        step_vals = [step_vals;t];
        x = x+t*direction;
        grad_x = grad(x,A);
        hess_x_inv = pinv(hess(x,A));
        iter = iter+1;
        lambdaSquare = grad_x'*hess_x_inv*grad_x;
    end
    figure(1)
    plot(iter_vals,f_vals)
    figure(2)
    stem(iter_vals,step_vals)
    figure(3)
    semilogy(iter_vals,eta_vals)
    f(x,A)
end
function [] = gradDescent()
    n = 100;
    m = 200;
    randn('state',1);
    A=randn(m,n);
    alpha = .01;
    beta = .5;
    x = randn(n,1)*.001;
    f_vals = [];
    iter_vals = [];
    eta_vals = [];
    step_vals = [];
    while(~dom(x,A))
        x = randn(n,1)*.001;
    end
    eta = 0.00001;
    grad_x = grad(x,A);
    
    maxIter = 1000;
    iter = 1;
    
    while(grad_x'*grad_x > eta && iter<maxIter)
        direction = -grad_x;
        t = 1;
        iter_vals = [iter_vals;iter];
        f_vals = [f_vals;f(x,A)];
        eta_vals = [eta_vals;grad_x'*grad_x];
        
        while(~dom(x+t*direction,A))
                t = beta*t;
        end
        while(f(x+t*direction,A) > f(x,A) + alpha*t*grad_x'*direction)
            t = beta*t;
        end
        step_vals = [step_vals;t];
        x = x+t*direction;
        grad_x = grad(x,A);
        iter = iter+1;
    end
    figure(1)
    plot(iter_vals,f_vals)
    figure(2)
    stem(iter_vals,step_vals)
    figure(3)
    semilogy(iter_vals,eta_vals)
end
function [f_x] = f(x,A)
    f_x = -sum(log(1-(A*x))) - sum(log(1-x.*x));
end
function [grad_x] = grad(x,A)
    grad_x = sum(diag((1./(1-A*x)))*A)' + 2*x./(1-x.^2);
end
function [hess_x] = hess(x,A)
    hess_x = (diag((1./(1-A*x)))*A)'*(diag((1./(1-A*x)))*A) + diag(2*(1+x.^2)./(1-x.^2).^2);
    %d = 1./(1-A*x);
    %hess_x = A'*diag(d.^2)*A + diag(1./(1+x).^2 + 1./(1-x).^2);
end
function [condition] = dom(x,A)
    condition = 0;
    if(sum(abs(x)<1) == length(x))
        if(sum(A*x<=1) == size(A,1))
            condition = 1;
        end
    end
end

        