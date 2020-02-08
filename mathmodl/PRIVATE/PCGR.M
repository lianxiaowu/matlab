function[p,posdef,k] = pcgr(DM,DG,g,kmax,tol,mtxmpy,fdata,H,R,pR);
%PCGR	Preconditioned conjugate gradients
%
% [p,posdef,k] = PCGR(DM,DG,g,kmax,tol,mtxmpy,fdata,H,R,pR) apply
% a preconditioned conjugate gradient procedure to the quadratic
% 
%         q(p) = .5p'Mp + g'p, where
%
% M = DM*H*DM + DG. kmax is a bound on the number of permitted
% CG-iterations, tol is a stopping tolerance on the residual (default
% is tol = .1), mtxmpy is the function that computes products
% with the Hessian matrix H (possible using data fdata),  
% and R is the cholesky factor of the preconditioner (transpose) of
% M. So, R'R approximates M(pR,pR), where pR is a permutation vector.
% On ouput p is the computed direction, posdef = 1 implies
% only positive curvature (in M) has been detected; posdef = 0
% implies p is a direction of negative curvature (for M).
% Output parameter k is the number of CG-iterations used (which
% corresponds to the number of multiplications with H).
%

%   Copyright (c) 1990-98 by The MathWorks, Inc.
%   $Revision: 1.3 $  $Date: 1998/08/17 19:42:07 $

% Initializations.
n = length(DG);
r = -g; 
p = zeros(n,1); 
val = 0;
m = 0;

% Precondition .
z = preproj(r,R,pR);
znrm = norm(z); 
stoptol = tol*znrm;
inner2 = 0; 
inner1 = r'*z;
posdef = 1;

kmax = max(kmax,1);  % kmax must be at least 1
% PRIMARY LOOP.
for k = 1:kmax
   if k==1
      d = z;
   else
      beta = inner1/inner2;
      d = z + beta*d;
   end
   ww = DM*d; 
   w = feval(mtxmpy,ww,H,fdata);
   ww = DM*w +DG*d;
   denom = d'*ww;
   if denom <= 0
      p = d/norm(d);
      posdef = 0;
      break
   else
      alpha = inner1/denom;
      p = p+ alpha*d;
      r = r - alpha*ww;
   end
   z = preproj(r,R,pR);
   
   % Exit?
   if norm(z) <= stoptol 
      break; 
   end
   inner2 = inner1;
   inner1 = r'*z;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function w = preproj(r,RPCMTX,ppvec);
%PREPROJ Apply preconditioner
%
% w = preproj(r,RPCMTX,ppvec) Apply a preconditioner to vector r.
% The conceptual preconditioner is H(ppvec,ppvec) = RPCMTX'*RPCMTX


% Initialization
n = length(r);
one = ones(n,1); 
vtol = 100*eps*one; 

if nargin < 3 | isempty(ppvec)
   ppvec = (1:n); 
   if nargin < 2 | isempty(RPCMTX) 
      RPCMTX = speye(n); 
   end
end 

%  Precondition
wbar = RPCMTX'\r(ppvec);
w(ppvec,1) = RPCMTX\wbar;




