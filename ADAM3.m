function [M,fiit] = ADAM3(prec, X,nc,varargin)
nd = ndims(X);
sz = size(X);
tsz = prod(sz);

params = inputParser;
params.addParameter('maxiters',10,@isscalar);
params.addParameter('festtol',-Inf,@isscalar);
params.addParameter('maxfails',50,@isscalar);
params.addParameter('beta1',0.9,@isscalar);
params.addParameter('beta2',0.999,@isscalar);
params.addParameter('epsilon',1e-8,@isscalar);
params.addParameter('rate', 1e-3, @isscalar);
params.addParameter('epciters', 10, @isscalar);
params.addParameter('decay', 0.1, @isscalar);
params.parse(varargin{:});

epciters = params.Results.epciters;
rate = params.Results.rate;
decay = params.Results.decay;
festtol = params.Results.festtol;
maxiters = params.Results.maxiters;
maxfails = params.Results.maxfails;
beta1 = params.Results.beta1;
beta2 = params.Results.beta2;
epsilon = params.Results.epsilon;


vecsz = sum(sz)*nc;
[fh,gh,lb]=gcp_setup3(prec, X);

m = cell(nd,1);
v = cell(nd,1);
M0 = cell(nd,1);
for k = 1:nd
    m{k}=zeros(sz(k),nc);
    v{k}=zeros(sz(k),nc);
    M0{k}=rand(sz(k),nc);
end
M0 = ktensor(M0);
fest = fg3(prec, M0,X,fh);


M = M0;
nfails = 0;
titers = 0;
Msave = M;
msave = m;
vsave = v;
fest_prev = fest;

fest_trace = zeros(maxiters+1,1);
step_trace = zeros(maxiters+1,1);
fest_trace(1)=fest;


for nepoch = 1:maxiters
    step = decay^nfails*rate;
    for iter = 1:epciters
        titers = titers + 1;
        [~,Gest] = fg3(prec,M,X,fh,gh);
        m = cellfun(@(mk,gk) beta1*mk + (1-beta1)*gk,m,Gest,'UniformOutput',0);
        v = cellfun(@(vk,gk) beta2*vk + (1-beta2)*gk.^2,v,Gest,'UniformOutput',0);
        mhat = cellfun(@(mk) mk/(1-beta1^titers),m,'UniformOutput',0);
        vhat = cellfun(@(vk) vk/(1-beta2^titers),v,'UniformOutput',0);
        M.u = cellfun(@(uk,mhk,vhk) max(lb,uk-step*mhk./(sqrt(vhk)+epsilon)),M.u,mhat,vhat,'UniformOutput',0);
    end
    fest = fg3(prec, M,X,fh);
    fest_trace(nepoch+1) = fest;
    step_trace(nepoch+1)=step;
    failed_epoch = fest>fest_prev;
    if failed_epoch
        nfails = nfails+1;
    end
    festtol_test = fest<festtol;
    if failed_epoch
        
        % Back up to best solution so far!
        M = Msave;
        m = msave;
        v = vsave;
        fest = fest_prev;
        titers = titers - epciters;
    else
        
        % Save current solution
        Msave = M;
        msave = m;
        vsave = v;
        fest_prev = fest;
        
    end
    if (nfails > maxfails) || festtol_test
        break;
    end
end
P = ktensor(M);
iprod = innerprod(X,P);
normresidual = sqrt( norm(X)*norm(X) + norm(P)*norm(P) - 2 * iprod );
fiit = 1 - (normresidual / norm(X));
