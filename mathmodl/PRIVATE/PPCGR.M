function[s,posdef,k,znrm,ex,lambda,PT,LZ,UZ,pcolZ,PZ] = ppcgr(grad,H,Aeq,b,...
   options,verb,fdata,computeLambda,tolA,PT,...
   LZ,UZ,pcolZ,PZ);
%PPCGR	Preconditioned conjugate gradients with linear equalities.
%
% [s,posdef,k] = PPCGR(grad,H,A,b,options,fdata)  apply
% a preconditioned conjugate gradient procedure to the linearly
% constrained quadratic
%
%         q(s) = .5s'Hs + grad's: As = b.
%
% On output s is the computed direction, posdef = 1 implies
% only positive curvature (in H constrained to null(A))
% has been detected; posdef = 0
% implies s is a direction of negative curvature (for M in null(A)).
% Output parameter k is the number of CG-iterations used (which
% corresponds to the number of multiplications with H).
% A preconditioner will be based on the form
%
%            HH  AA'
%    M = 
%            AA  0
%
% where in general HH is a SPD banded approximation
% to H, the nonzeros of AA are a subset of the nonzeros of A (based
% on a stopping tolerance). Parameters in the named parameter
% list options can be used to overide default values to define M,
% including defining M (H) implicitly.
%
% [s,posdef,k] = PPCGR(grad,H,A,b,options,fdata,tolA) 
% Overide the dropping tolerance to define AA from A.
%
% [s,posdef,k] = PPCGR(grad,H,A,b,options,fdata,tolA,PT,LZ,UZ,pcolZ,PZ);
% PT is a permutation matrix such that A*PT has a leading nonsingular 
% m-by-m matrix, A_1, where m is the number of rows of A
% PZ*A_1(:,pcolZ) = LZ*UZ, the sparse LU-factorization on A_1.
%
% [s,posdef,k,znrm,ex,PT,LZ,UZ,pcolZ,PZ]=PPCGR(grad,H,A,b,options,fdata) 
% ouputs the sparse LU-factorization of A_1. znnrm is the scaled norm
% of the residual, ex is the termination code.

%   Copyright (c) 1990-98 by The MathWorks, Inc.
%   $Revision: 1.12 $  $Date: 1998/09/09 19:48:00 $


if nargin < 3, 
   error('ppcgr requires at least 3 input arguments'); 
end
[m,n] = size(Aeq); 
ex = 1;
s = [];
posdef = [];
k = [];
znrm = [];

if m > n, 
   ex = -2; s = []; lambda = []; posdef=[];k=[];znrm=[];
   %error(' ppcgr requires m < n');
   return
end

if nargin < 4 | isempty(b), 
   b = zeros(m,1); 
end
[mm,nn] = size(b);

if mm ~= m, 
   error('The number of rows in Aeq is incompatible with the length of beq.'); 
end

if nargin < 8,
   computeLambda = 0;
   if nargin < 7,
      fdata = []; 
      if nargin < 6, 
         verb = 0;
         if nargin < 5
            options = [];
         end, end, end, end

% Later Preconditioner and HessMult will be user-settable
pcmtx = optimget(options,'Preconditioner','preaug'); 
mtx = optimget(options,'HessMult','hmult') ;

pcf = optimget(options,'PrecondBandWidth') ;
if nargin >= 9 & ~isempty(tolA), 
   pcf(2,1) = tolA; % tolA not currently used 
end
tol = optimget(options,'TolPCG', 0.1) ;
kmax = optimget(options,'MaxPCGIter', n) ;
% In case the defaults were gathered from calling: optimset('quadprog'):
numberOfVariables = n;
if ischar(kmax)
   kmax = eval(kmax);
end
kmax = max(kmax,1);   % kmax must be at least 1

% save initial grad for computing Lagrange multipliers
gradsave = grad;

% Get feasible point
y = zeros(n,1);
if norm(b) > 1e-12
   [y, flag] = feasibl(Aeq,b);
   % We don't yet handle this case of dependent rows
   if flag == -1
      ex = -2; s = []; lambda = []; posdef=[];k=[];znrm=[];
      return
   end
   
   w = feval(mtx,y,H,fdata);
   grad = grad + w;
end

% Transform to get fundamental null basis
if nargin < 10
   PT = findp(Aeq);
end

A = Aeq*PT'; 
grad = PT*grad;

% Project grad
if nargin < 14 % last 4 arguments are a group
   [g,LZ,UZ,pcolZ,PZ] = fzmult(A,grad,'transpose');
else
   g = fzmult(A,grad,'transpose',LZ,UZ,pcolZ,PZ);
end
r = -g; 

nz = length(g); 
s = zeros(nz,1); 
val = 0;

% Compute preconditioner
HH = [];
if ~isempty(H);
   HH = PT*H*PT';
end
[L,U,P,pcol] = feval(pcmtx,HH,pcf,A);

% Precondition .
rr = [zeros(m,1);r];
rhs = [rr;zeros(m,1)];
zz = L\(P*rhs); 
zt(pcol,1) = U\zz; 
z = zt(m+1:n,1);
znrm = norm(z); 
stoptol = tol*znrm;
inner2 = 0; 
inner1 = r'*z; 
posdef = 1;

% PRIMARY LOOP.
for k = 1:kmax
   if k==1
      d = z;
   else
      beta = inner1/inner2; 
      d = z + beta*d;
   end
   dd = fzmult(A,d,'',LZ,UZ,pcolZ,PZ);
   dd = PT'*dd;
   w = feval(mtx,dd,H,fdata);
   ww = PT*w; 
   w = fzmult(A,ww,'transpose',LZ,UZ,pcolZ,PZ);
   denom = d'*w;
   if denom <= 0
      if norm(d) > 0
         s = d/norm(d);
      else
         s = d;
      end
      s = y + PT'*fzmult(A,s,'',LZ,UZ,pcolZ,PZ);
      posdef = 0; 
      if computeLambda
         w = feval(mtx,s,H,fdata); 
         rhs = -gradsave-w;
         lambda.eqlin = Aeq'\rhs;
         lambda.ineqlin = []; lambda.lower = []; lambda.upper = [];
      else
         lambda = [];
      end
      if denom == 0 & sum(d)== 0 % No curvature left: we have a (singular) solution
         ex = 1;
         if verb > 0
            disp('Optimization terminated successfully:')
            disp(' Local minimum found; the solution is singular.');
         end
      else % either denom < 0, or denom = 0 and sum(d) ~= 0   
         ex = -1; % this was ex = 3;
         if verb > 0
            if denom < 0
               disp('Exiting: Negative curvature direction detected. ')
               disp('         The solution is unbounded and at infinity;')
               disp('         Constraints are not restrictive enough.');
            else
               disp('Exiting:')
               disp('         Zero curvature direction detected.');
               disp('         The solution is unbounded and at infinity;')
               disp('         Constraints are not restrictive enough.');
            end
         end % verb > 0
      end
      return
   else % denom > 0, continue
      alpha = inner1/denom;
      s = s+ alpha*d; 
      r = r - alpha*w;
   end
   
   % Precondition
   rr =[zeros(m,1);r];  
   rhs = [rr;zeros(m,1)];
   zz = L\(P*rhs); 
   zt(pcol,1) = U\zz; 
   z = zt(m+1:n,1);
   
   % Exit?
   znrm = norm(z);
   if znrm <= stoptol,
      s = y +  PT'*fzmult(A,s,'',LZ,UZ,pcolZ,PZ); 
      if computeLambda
         w = feval(mtx,s,H,fdata); 
         rhs = -gradsave-w;
         lambda.eqlin = Aeq'\rhs;
         lambda.ineqlin = []; lambda.lower = []; lambda.upper = [];
      else 
         lambda = [];
      end   
      ex = 1;
      if verb > 0
         disp('Optimization terminated successfully:')
         disp(' Relative (projected) residual of PCG iteration <= OPTIONS.TolPCG');
      end
      return;
   end
   inner2 = inner1; 
   inner1 = r'*z;
end

if k >= kmax, 
   ex = 0;  % zero is what it should be, was ex = 2;
   if verb > 0
      disp('Maximum number of PCG iterations exceeded;')
      disp('   increase options.MaxPCGIter')
   end              
end
s = y +  PT'*fzmult(A,s,'',LZ,UZ,pcolZ,PZ);
if computeLambda
   w = feval(mtx,s,H,fdata);  
   rhs = -gradsave-w;
   lambda = Aeq'\rhs;
   lambda.ineqlin = []; lambda.lower = []; lambda.upper = [];
else 
   lambda = [];
end


