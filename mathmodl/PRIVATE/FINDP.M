function P  = findp(A);
%FINDP  Nonsingular basis permutation matrix.
%   P = FINDP(A) locates a permutation matrix such that
%   the leading square matrix is nonsingular.If A has 
%   more columns (n) than rows(m), and A is of rank m, then
%   P is a permutation matrix such that the leading m-matrix
%   of A*P is of full rank.

%   Copyright (c) 1990-98 by The MathWorks, Inc.
%   $Revision: 1.2 $  $Date: 1998/03/21 16:29:08 $

[m,n] = size(A);
if m < n, 
   A = A';
end
if ~issparse(A), 
   A=sparse(A); 
end;
pp = colmmd(A);
[L,U,P] = lu(A(:,pp));


