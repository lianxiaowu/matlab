function A = gangstr(M,tol);
%GANGSTR Replace small nonzeroes with zero.
%   A = GANGSTR(M,tol) replaces nonzero entries of M that
%   are less than tol in magnitude with zero values, subject to
%   maintenance of full structural rank. tol is increased by
%   factors of 10, if needed, until strank(A) = strank(M).
%
%   GANGSTR expects A to have more columns than rows.

%   Copyright (c) 1990-98 by The MathWorks, Inc.
%   $Revision: 1.2 $  $Date: 1998/03/27 01:29:29 $

if isempty(M)
   A=M;
   return
end

[m,n] = size(M);
if m > n
   error('M matrix must have more columns than rows.')
end

if nargin < 2, 
   tol = 1e-2; 
end

% Normalize M
Msqr=M.*M;
% Sum rows of Msqr and return column vector.
% Msqr cannot be a column vector.
X=(sum(Msqr'))'; % spdiags needs a column vector as input

X(X==0) = 1;
M=spdiags(1./sqrt(X),0,m,m)*M;

% Remove nonzeros
dim = sprank(M);
sprA = 0;
A=M;
[I,J,V]=find(M);
while sprA < dim
   % Take the max of the columns of abs(M) unless M is one row only   
   if m > 1 
      maxvec=full(tol*max(abs(M))); 
      maxvec=maxvec(:);    % make maxvec consistent with I,J,V
   else 
      maxvec=full(tol*abs(M)); % maxvec, I,J,V are all rows in this case 
   end
   tobekept=find(abs(V) > maxvec(J));
   A=sparse(I(tobekept),J(tobekept),V(tobekept),m,n);
   sprA = sprank(A);
   tol = tol/10;
end

% Denormalize
A=spdiags(sqrt(X),0,m,m)*A;
