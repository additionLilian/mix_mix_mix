function [u,e,T] = SGD_outer3(prec, X,iterCG,iterSG,R,step)
s = size(X);
N = ndims(X);
A = cell(N,1);
e = [];
rng(15);
for i = 1:N
     A{i} = rand(s(i),R);
end
tic;
T = zeros(1,iterSG);
for i = 1:iterSG
    [A,error] = CG_cp3(prec,X,A,iterCG,step);
    e(end+1)=error;
    T(i)=toc;
end
u = A;
normX = norm(X);
e = e/normX;