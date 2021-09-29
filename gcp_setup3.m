function [fh,gh,lowerbnd] = gcp_setup3(prec,X)

% suppose: normal distribution 
if prec == 0
    fh = @(x,m) (half(m)-half(x)).^2;
    gh = @(x,m) 2.*(half(m)-half(x));
elseif prec == 1
    fh = @(x,m) (single(m)-single(x)).^2;
    gh = @(x,m) 2.*(single(m)-single(x));
elseif prec == 2
    fh = @(x,m) (m-x).^2;
    gh = @(x,m) 2.*(m-x);
end
lowerbnd = -Inf;

function tf = valid_nonneg(X)

if isa(X,'sptensor')
    tf = all(X.vals > 0);
else
    tf = all(X(:) > 0);
end

function tf = valid_binary(X)

if isa(X,'sptensor')
    tf = all(X.vals == 1);
else
    tf = isequal(unique(X(:)),[0;1]);
end

function tf = valid_natural(X)

if isa(X, 'sptensor')
    vals = X.vals;
else
    vals = X(:);
end

tf = all(vals >= 0) && all(vals == round(vals));
