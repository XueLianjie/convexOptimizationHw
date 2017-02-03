L = 0;
U = 1;
k = 201;
t = [-3:abs(6)/200:3]';
y = exp(t);
%x = [a0;a1;a2;b1;b2]
C = [ones(k,1) t t.^2 -y.*t -y.*t.^2];
F = [zeros(k,1) zeros(k,1) zeros(k,1) t t.^2];
g = ones(k,1);
d = -y;
u=U;
l=L;
for j = 1:ceil(log2((U-L)/.001))
    alpha = (u+l)/2;
    cvx_begin
        variables x(5)
        subject to
            abs(C*x + d) <= alpha*(F*x + g)
            F*x + g >= 0
    cvx_end
    if strcmp(cvx_status,'Solved')
        u = alpha;
        x_opt = x;
        alpha_opt = alpha;
    else
        l = alpha;
    end
end
a0 = x_opt(1);
a1 = x_opt(2);
a2 = x_opt(3);
b1 = x_opt(4);
b2 = x_opt(5);

y_fit = (a0 + a1*t + a2*t.^2) / (1 + b1*t + b2*t.^2);
figure(1)
plot(t,y);
figure(2)
plot(t,y_fit);
        
        