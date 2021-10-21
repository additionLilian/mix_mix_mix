function [A,error] = CG_cp3(prec,X,initial,maxiter,step)
rng(8);
M = 40;
N = ndims(X);
A = initial;
U = ktensor(A);
error = [norm(minus(full(U),full(X)))];
G = gradient_calc2(2,X,A);
D = cellfun(@(x)(-x),G,'UniformOutput',0);
s = size(X);
for i = 1:maxiter
    step = max(step/1.5,0.1);
    ind = [];
    for j = 1:N
        tmp = randi(s(j),M,1);
        ind = [ind,tmp];
    end
    Aold = A;
    Dold = D;
    U = ktensor(A);
    error_old = norm(minus(full(U),full(X)));
    Dold = cellfun(@(x)(double(x)),Dold,'UniformOutput',0);
    Dnorm = sum(cellfun(@(x)norm(x), Dold));
    Dteld = cellfun(@(x)(x/Dnorm),Dold,'UniformOutput',0);
    AsD = cellfun(@(x,y)x+y*(step),A,Dteld,'UniformOutput',0);
    GsD = gradient_M(prec,AsD,X,ind);
    %GsD = gradient_calc2(prec, tensor(X),AsD,R,ind);
    Gold_cell = cell(1,N);
    Gnew_cell = cell(1,N);
    hessian_cell = cell(1,N);
    vecD = cell(1,N);
    for k = 1:N
        [a,b] = size(G{k});
        Gold_cell{k} = double(reshape(G{k},a*b,1));
        vecAk = reshape(A{k},a*b,1);
        vecD{k} = double(reshape(D{k},a*b,1));
        vecGsD = double(reshape(GsD{k},a*b,1));
        hessian_cell{k} = ((vecGsD - Gold_cell{k})/step)*Dnorm;
        alpha = -Gold_cell{k}'*vecD{k}/(vecD{k}'*hessian_cell{k});
        vecAk = vecAk + alpha * vecD{k};
        Ak = reshape(vecAk,a,b);
        A{k} = Ak;
    end
    G = gradient_M(prec,A,X,ind);
    %G = gradient_calc2(prec,X,A,R,ind);
    for k = 1:N
        [a,b] = size(G{k});
        Gnew_cell{k} = reshape(G{k},a*b,1);
        beta = norm(double(Gnew_cell{k}),2)^2/norm3(double(Gold_cell{k}))^2;
        
        vecD{k} = -Gnew_cell{k}+beta*vecD{k};
    end
    A = cellfun(@(x)(double(x)),A,'UniformOutput',0);
    U = ktensor(A);
    
    error_new = norm(minus(full(U),full(X)));
    if error_new>error_old
        A = Aold;
        D = Dold;
    else
        Gold_cell = Gnew_cell;
        error(end+1) = error_new;
    end
end
error = error(end);