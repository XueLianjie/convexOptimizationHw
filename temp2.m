mu1 = 8;
mu2 = 20;
sigma1 = 6;
sigma2 = 17.5;
rho = -.25;
n = 100;
r = [-30:(100/(n-1)):70]';
p1 = exp(-(r-mu1).^2./(2*sigma1.^2))./sum(exp(-(r-mu1).^2./(2*sigma1.^2)));
p2 = exp(-(r-mu2).^2./(2*sigma2.^2))./sum(exp(-(r-mu2).^2./(2*sigma2.^2)));
r_sum = repmat(r,1,100);
r_sum = r_sum + r_sum';
r_prod = r_sum.*r_sum';
r_mask = r_sum<=0;

cvx_begin
    variables P(n,n)
    
    minimize (ones(n,1)'*(r_mask.*P)*ones(n,1))
    subject to
        P*ones(n,1) == p2;
        P'*ones(n,1) == p1;
        (r-mu2)'*P*(r-mu1) == rho*sigma1*sigma2;
        P >= 0;
        ones(n,1)'*P*ones(n,1) == 1
cvx_end

mesh(P)
