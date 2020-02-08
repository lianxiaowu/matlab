% $Revision: 1.2 $
function[L,U,P,pcol] = preaug(H,pcf,A);
%PREAUG  example preconditioner of augmented matrix
%
% [L,U,P,pcol] = PREAUG(H,pcf,A) computes a sparse
% factorization of the LU-factorization of
%
%                   H     A'
%        M =          
%                   A     0  
%
%
% say,
%
%                   HH   AA'
%        MM =   
%                   AA    0
%
% where HH is SPD, usually banded, the nonzeros of AA
% are a subset of the nonzeros of A. If 0 < pcf(1,1) < n then
% the upper bandwidth of HH is pcf(1,1); if pcf(1,1) >= n then HH
% is a sparse Cholesky factorization of H. pcf(2,1) is the dropping
% tolerance for A --> AA.
%
%

%   Copyright (c) 1990-98 by The MathWorks, Inc.
%   $Revision: 1.2 $  $Date: 1998/03/21 16:29:03 $
 
%
% Approximate H
if nargin ==0, 
   error('no input parameters in function preaug'); 
end
if isempty(H), 
   error('input parameter H cant be null in function preaug'); 
end
n = length(H);
epsi = .001*ones(n,1);
info = 1;
HH = H;
if nargin <2, 
   pcf = [0;.05]; 
end
[mpcf,dumm] = size(pcf);
if mpcf < 2 
   pcf(2,1) = .05; 
end
tolA = abs(pcf(2,1));
ppcf = pcf(1,1);

if ppcf >= n % Try complet approx to H
   p = symmmd(H);
   ddiag = diag(H);
   mind = min(ddiag);
   lambda = 0;
   if mind < 0, 
      lambda = -mind + .001; 
   end
   while info > 0
      H = H + lambda*speye(n);
      [R,info] = chol(H(p,p));
      lambda = lambda + 10;
   end
elseif (ppcf > 0) & ( ppcf < n) %Banded approximation to H
   % Cluster diagonal
   p = symrcm(H);
   HH = H(p,p);
   bndw = ppcf;
   HH = tril(triu(HH,-bndw),bndw);
   lambda = 0;
   ddiag = diag(HH);
   mind = min(ddiag);
   if mind < 0, 
      lambda = -mind + .001; 
   end
   while info > 0
      HH = HH + lambda*speye(n);
      [R,info] = chol(HH);
      lambda = 4*lambda;
      if lambda <= .001, 
         lambda = 1; 
      end
   end
   H(p,p) = HH;
else % diagonal approximation for H
   dnrms = sqrt(sum(H.*H))';
   d = max(dnrms,epsi);
   H = sparse(1:n,1:n,full(d));
end
if nargin < 3, 
   A = []; 
end
if isempty(A)
   L = R'; 
   U = R; 
   pcol = p; 
   P = sparse(n,n,pcol);
   return
end
[m,nn] = size(A);
if nn ~= n, 
   error('dim. of A conflict with d. of H in fcn preaug'); 
end
if m > n, 
   error('matrix A should have more cols than rows'); 
end
%
% Approximate A
%
% Set up the augmented matrix
A = gangstr(A,tolA);
M = [H,A';A, sparse(m,m)];
pcol = colmmd(M); 
pcol = pcol;
[L,U,P] = lu(M(:,pcol));


