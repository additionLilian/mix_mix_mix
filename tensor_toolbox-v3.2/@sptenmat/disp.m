function disp(a,name)
%DISP Command window display of a sptenmat.
%
%   DISP(T) displays the tensor without printing its name.
%
%   DISP(T,NAME) displays the tensor with the given name.
%
%   See also SPTENMAT/DISPLAY.
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>



if ~exist('name','var')
    name = 'ans';
end

% Extract the number of nonzeros and number of dimensions
nz = size(a.vals,1);

% Print an intro sentence giving the name and the size
if (nz == 0)
    fprintf('%s is an all-zero sptenmat from an sptensor of size %s\n',...
        name, tt_size2str(a.tsize));
else
    fprintf('%s is a sptenmat from an sptensor of size %s with %d nonzeros\n',...
        name, tt_size2str(a.tsize), nz);
end

fprintf(1,'\t%s.rindices = %s (modes of tensor corresponding to rows)\n',...
        name,['[ ' num2str(a.rdims) ' ]'] );
fprintf(1,'\t%s.cindices = %s (modes of tensor corresponding to columns)\n',...
        name,['[ ' num2str(a.cdims) ' ]'] );


% Stop insane printouts
if (nz > 1000)
    r = input('Are you sure you want to print all nonzeros? (Y/N) ','s');
    if upper(r) ~= 'Y', return, end;
end

% Return now if there are no nonzeros
if (nz == 0)
    return;
end

% Pre-allocate memory for the output
output = cell(nz,1);
spc = floor(log10(max(a.subs)))+1;
fmt = ['\t(%' num2str(spc(1)) 'd,%' num2str(spc(2)) 'd)\t%g'];

for i = 1:nz
    output{i} = sprintf(fmt, a.subs(i,:), a.vals(i));
end
fprintf('%s\n',output{:});



