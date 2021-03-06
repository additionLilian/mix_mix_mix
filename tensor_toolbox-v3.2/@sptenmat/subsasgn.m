function t = subsasgn(t,s,b)
%SUBSASGN Subscripted assignment for sptenmat.  
%
%   Examples 
%   X = sptenmat(sptenrand([3 4 2],10),1);
%   X(1:2,1:2) = ones(2,2); <-- Calls SUBSASGN 
%
%   See also SPTENMAT, SPTENMAT/SUBSREF.
%
%Tensor Toolbox for MATLAB: <a href="https://www.tensortoolbox.org">www.tensortoolbox.org</a>



%TODO: Make implementation efficient. It's not right now.
%TODO: Add error checking.

switch s.type    
    case '()'
        rsubs = s.subs{1};
        csubs = s.subs{2};

        % If called with single element in rhs, then create a vector
        if numel(b) == 1
          b = repmat(b, numel(rsubs) * numel(csubs), 1);
        end

        % Initialize some variables for new entries in the matrix
        newsubs = [];
        newvals = [];

        k = 0;

        % Loop over the row and column indices, finding the
        % appropriate row index for the (i,j) subscript
        for j = 1:length(csubs)
          indxc = find(t.subs(:,2) == csubs(j));
          for i = 1:length(rsubs)
            indxr = find(t.subs(indxc,1) == rsubs(i));
            indx = indxc(indxr);

            k = k + 1;  % increment counter into b
            if isempty(indx)
              newsubs = [newsubs; rsubs(i) csubs(j)];
              newvals = [newvals; b(k)];
            else
              %t.subs(indx,:);
              t.vals(indx) = b(k);
            end
          end
        end

        % If there are new values to append, then add them on and sort
        if ~isempty(newvals)
          t.subs = [t.subs; newsubs];
          t.vals = [t.vals; newvals];
          [t.subs,indx] = sortrows(t.subs);
          t.vals = t.vals(indx);
        end

    otherwise
        error('Invalid assignment for sptenmat.')
end
