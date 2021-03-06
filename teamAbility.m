train=[1 2 1;
1 3 1;
1 4 1;
1 5 1;
1 6 1;
1 7 1;
1 8 1;
1 9 1;
1 10 1;
2 3 -1;
2 4 -1;
2 5 -1;
2 6 -1;
2 7 -1;
2 8 -1;
2 9 -1;
2 10 -1;
3 4 1;
3 5 -1;
3 6 -1;
3 7 1;
3 8 1;
3 9 1;
3 10 1;
4 5 -1;
4 6 -1;
4 7 1;
4 8 1;
4 9 -1;
4 10 -1;
5 6 1;
5 7 1;
5 8 1;
5 9 -1;
5 10 1;
6 7 1;
6 8 1;
6 9 -1;
6 10 -1;
7 8 1;
7 9 1;
7 10 -1;
8 9 -1;
8 10 -1;
9 10 1;
];

test=[1 2 1;
1 3 1;
1 4 1;
1 5 1;
1 6 1;
1 7 1;
1 8 1;
1 9 1;
1 10 1;
2 3 -1;
2 4 1;
2 5 -1;
2 6 -1;
2 7 -1;
2 8 1;
2 9 -1;
2 10 -1;
3 4 1;
3 5 -1;
3 6 1;
3 7 1;
3 8 1;
3 9 -1;
3 10 1;
4 5 -1;
4 6 -1;
4 7 -1;
4 8 1;
4 9 -1;
4 10 -1;
5 6 -1;
5 7 1;
5 8 1;
5 9 1;
5 10 1;
6 7 1;
6 8 1;
6 9 1;
6 10 1;
7 8 1;
7 9 -1;
7 10 1;
8 9 -1;
8 10 -1;
9 10 1;
];

n = 10;

m = 45;

m_test = 45;

sigma= 0.250;

A = zeros(m,n);
for i=1:m
    for j=1:m
        if(j == train(i,1))
            A(i,j) = train(i,3);
        elseif(j == train(i,2))
            A(i,j) = -train(i,3);
        end
    end
end
y_train = train(:,3);

cvx_begin
    variables a(n)
    v = (1/sigma).*A*a;
    maximize(sum(log_normcdf(v)))
    subject to
        a >= 0
        a <= 1
cvx_end

y_test = test(:,3);
y_pred = sign(a(test(:,1))-a(test(:,2)));

sum(y_test == y_pred)/length(y_pred)
sum(y_test == y_train)/length(y_pred)