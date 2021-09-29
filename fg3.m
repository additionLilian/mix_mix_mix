function [F,G] = fg3(prec, M, X, f, g, W, computeF, computeG, vectorG)
Mfull = full(M);
Mv = Mfull(:);
Xfull = full(X);
Xv = Xfull(:);
F = [];
G = [];

if nargin < 8  
    
    if ~exist('W','var')
        W = [];
    end
    
    if ~exist('computeF','var')
        computeF = true;
    end
    
    if ~exist('computeG','var') 
        computeG = (nargout > 1);
    end
    
    if ~exist('vectorG','var')
        vectorG = false;
    end
    
end

%% Calculate function value
if computeF
    
    Fvec = f(Xv, Mv); % F is a vector
    
    if ~isempty(W)
        Fvec = W(:).*Fvec;  % be sure to zero out any unknown entries
    end
    
    F = sum(Fvec);
    
end

%% QUIT IF ONLY NEED FUNCTION EVAL
if ~computeG
    return;
end

%% Gradient calculation
Y = g(Xv,Mv); % Result is a vector
Y = tensor(Y,size(X));
if ~isempty(W)
    Y = W.*Y;
end
%% Gradient wrt U's using MTTKRP sequence.
G = mttkrps3(prec,Y,M.u);
%% Assemble gradient 
if vectorG
    G = cell2mat(cellfun(@(x) x(:), G, 'UniformOutput', false));
end

%% If not computing F, set F (the 1st return arugment) to be the gradient
if ~computeF
    F = G;
end
