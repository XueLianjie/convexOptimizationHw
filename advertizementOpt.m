% data for online ad display problem
rand('state',0);
n=100;      %number of ads
m=30;       %number of contracts
T=60;       %number of periods

I=10*rand(T,1);  %number of impressions in each period
R=rand(n,T);    %revenue rate for each period and ad
q=T/n*50*rand(m,1);     %contract target number of impressions
p=rand(m,1);  %penalty rate for shortfall
Tcontr=(rand(T,m)>.8); %one column per contract. 1's at the periods to be displayed
for i=1:n
	contract=ceil(m*rand);
	Acontr(i,contract)=1; %one column per contract. 1's at the ads to be displayed
end

%{
cvx_begin
    variables N(n,T)
    totalRev = sum(sum(N.*R));
    totalPenalty = 0;
    for i=1:m
        totalPenalty = totalPenalty + p(i)*max(q(i)-sum(sum(N.*(Acontr(:,i)*Tcontr(:,i)'))),0)
    end
    maximize (totalRev-totalPenalty)
    subject to
        N >= 0
        sum(N)' == I
cvx_end
%}
%if N is chosen so that ads with highest revenue are shown
N = [];
for i=1:n
    N = [N;max(R)];
end
N = (N==R)*diag(I);
totalRev = sum(sum(N.*R));
totalPenalty = 0;
for i=1:m
    totalPenalty = totalPenalty + p(i)*max(q(i)-sum(sum(N.*(Acontr(:,i)*Tcontr(:,i)'))),0)
end
totalRev
totalPenalty
profit = totalRev-totalPenalty
    