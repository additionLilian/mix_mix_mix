function P = mykhatrirao3(prec,A,order)

if order == 'r'
    matorder = length(A):-1:1;
else
    matorder = 1:length(A);
end

%% Error check on matrices and compute number of rows in result 
ndimsA = cellfun(@ndims, A);
ncols = cellfun(@(x) size(x, 2), A);
if(~all(ncols == ncols(1)))
    error('All matrices must have the same number of columns.');
end


%% Computation
if prec == 1
    A = cellfun(@(x)(single(x)),A,'UniformOutput',0);
elseif prec == 0
    A = cellfun(@(x)(half(x)),A,'UniformOutput',0);
end
N = ncols(1);
P = A{matorder(1)};
for i = matorder(2:end)
    [ma,na] = size(P);
    Ai = A{i};
    [mb,~] = size(Ai);
    if prec == 0
        PP = half(zeros(ma*mb,na));
        A{2} = A{2}./60000;
    elseif prec == 1
        PP = single(zeros(ma*mb,na));
    else
        PP = zeros(ma*mb,na);
    end
    for j = (1:na)
        PP(:,j) = double(kron(P(:,j),Ai(:,j)));
    end
    P = PP;
end
P = reshape(P,[],N);
