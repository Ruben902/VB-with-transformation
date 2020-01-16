function M = DuplicationM(n, option)
% function M = DuplicationM(n)
% M = DuplicationM(n, 'lo') % (default) OR
% M = DuplicationM(n, 'up') % 
% Return duplication matrix order n
%
% It is always assumed Duplication arount main diagonal (k=0)
%
% Ouput are sparse
%
% DuplicationM(size(A),'lo')*A(itril(size(A))) == A(:) %true for lower half
% DuplicationM(size(A),'up')*A(itriu(size(A))) == A(:) %true for upper half
%
% Author: Bruno Luong <brunoluong@yahoo.com>
% Date: 21/March/2009
%
% Ref : Magnus, Jan R.; Neudecker, Heinz (1980), "The elimination matrix:
% some lemmas and applications", Society for Industrial and Applied Mathematics.
% Journal on Algebraic and Discrete Methods 1 (4): 422449,  
% doi:10.1137/0601049, ISSN 0196-5212.

if nargin<2
    option = 'lo'; % default
end

if isscalar(n)
    n = [n n];
end

switch lower(option(1))
    case 'l' % u, lo, LO, LOWER ...
        [I,J] = itril(n);
    case 'u' % u, up, UP, UPPER ...
        [I,J] = itriu(n);
    otherwise
        error('option must be ''lo'' or ''up''');
end

% Find the sub/sup diagonal part that can flip to other side
loctri = find(I~=J & J<=n(1) & I<=n(2));
% Indices of the flipped part
Itransposed = sub2ind(n, J(loctri), I(loctri));

% Convert to linear indice
I =  sub2ind(n, I, J);

% Result
M = sparse([I; Itransposed], ...
           [(1:length(I))'; loctri], 1, prod(n), length(I));

end


function [I J] = itril(sz, k)
% function [I J] = itril(sz, k) % OR
% I = itril(sz, k)
%
% Return the subindices [I J] (or linear indices I if single output call)
% in the purpose of extracting an lower triangular part of the matrix of
% the size SZ. Input k is optional shifting. For k=0, extract from the main
% diagonal. For k>0 -> above the diagonal, k<0 -> below the diagonal
% 
% This returnd same as [...] = find(tril(ones(sz),k))
% - Output is a column and sorted with respect to linear indice
% - No intermediate matrix is generated, that could be useful for large
%   size problem
% - Mathematically, A(itril(size(A)) is called (lower) "half-vectorization"
%   of A 
%
% Example:
%
% A = [ 7     5     4
%       4     2     3
%       9     1     9
%       3     5     7 ]
%
% I = itril(size(A))  % gives [1 2 3 4 6 7 8 11 12]'
% A(I)                % gives [7 4 9 3 2 1 5  9  7]' OR A(tril(A)>0)
%
% Author: Bruno Luong <brunoluong@yahoo.com>
% Date: 21/March/2009

if isscalar(sz)
    sz = [sz sz];
end
m=sz(1);
n=sz(2);

% Main diagonal by default
if nargin<2
    k=0;
end

nc = min(n,m+k); % number of columns of the triangular part
lo = max((1:nc).'-k,1); % lower row indice for each column
hi = m + zeros(nc,1); % upper row indice for each column

if isempty(lo)
    I = zeros(0,1);
    J = zeros(0,1);
else
    c=cumsum([0; hi-lo]+1); % cumsum of the length
    I = accumarray(c(1:end-1), (lo-[0; hi(1:end-1)]-1), ...
                   [c(end)-1 1]);
    I = cumsum(I+1); % row indice
    J = cumsum(accumarray(c,1));
    J = J(1:end-1); % column indice
end

if nargout<2
    % convert to linear indices
    I = sub2ind([m n], I, J);
end

end % itril


function [I J] = itriu(sz, k)
% function [I J] = itriu(sz) % OR
% I = itriu(sz) OR
% 
% Return the subindices [I J] (or linear indices I if single output call)
% in the purpose of extracting an upper triangular part of the matrix of
% the size SZ. Input k is optional shifting. For k=0, extract from the main
% diagonal. For k>0 -> above the diagonal, k<0 -> below the diagonal
%
% This returnd same as [...] = find(triu(ones(sz),k))
% - Output is a column and sorted with respect to linear indice
% - No intermediate matrix is generated, that could be useful for large
%   size problem
% - Mathematically, A(itriu(size(A)) is called (upper) "half-vectorization"
%   of A 
%
% Example:
%
% A = [ 7     5     4
%       4     2     3
%       9     1     9
%       3     5     7 ]
%
% I = itriu(size(A))  % gives [1 5 6 9 10 11]'
% A(I)                % gives [7 5 2 4  3  9]' OR A(triu(A)>0)
%
% Author: Bruno Luong <brunoluong@yahoo.com>
% Date: 21/March/2009

if isscalar(sz)
    sz = [sz sz];
end
m=sz(1);
n=sz(2);

% Main diagonal by default
if nargin<2
    k=0;
end

nc = n-max(k,0); % number of columns of the triangular part
lo = ones(nc,1); % lower row indice for each column
hi = min((1:nc).'-min(k,0),m); % upper row indice for each column

if isempty(lo)
    I = zeros(0,1);
    J = zeros(0,1);
else
    c=cumsum([0; hi-lo]+1); % cumsum of the length
    I = accumarray(c(1:end-1), (lo-[0; hi(1:end-1)]-1), ...
                   [c(end)-1 1]);
    I = cumsum(I+1); % row indice
    J = accumarray(c,1);
    J(1) = 1 + max(k,0); % The row indices starts from this value
    J = cumsum(J(1:end-1)); % column indice
end

if nargout<2
    % convert to linear indices
    I = sub2ind([m n], I, J);
end

end % itriu




