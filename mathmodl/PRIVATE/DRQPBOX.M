function[s,posdef,pcgit] = drqpbox(D,DS,grad,delta,g,dv,mtxmpy,...
   data,pcmtx,pcoptions,tol,H,llsprob,kmax);
%DRQPBOX Descent direction for quadratic problem.
%	[s,posdef,pcgit] = DRQPBOX(D,DS,grad,delta,g,dv,mtxmpy,...
%                                    data,pcmtx,pcoptions,tol,H,llsprob)
%   determines s, a descent direction (for use with SLLSBOX,SQPBOX) 
%   for quadratic minimization subject to box constraints.
%   If negative curvature is discovered in the CG process
%   then posdef = 0; otherwise posdef = 1. pcgit is the
%   number of CG iterations. LLSPROB is a flag that tells if
%   the caller is SLLSBOX or SQPBOX so the preconditioner can be
%   called correctly.

%   Copyright (c) 1990-98 by The MathWorks, Inc.
%   $Revision: 1.3 $  $Date: 1998/08/17 19:42:06 $

n = length(grad); 
tsize = 0; 
pcgit = 0; 
tol2 = sqrt(eps);

DM = D*DS;
tau = 1e-4; 
vtau = tau*ones(n,1); 
diagDM = full(diag(DM));
ddiag = abs(g).*dv; 
arg = (abs(g) + diagDM < tau);
ddiag(arg == 1) = vtau(arg == 1);
DG = sparse(1:n,1:n,full(ddiag));

% A PRECONDITIONED CONJUGATE GRADIENT ROUTINE IS USED.
if llsprob
   [R,permR] = feval(pcmtx,DM,DG,H,pcoptions);
else
   HH = DM*H*DM + DG;
   [R,permR] = feval(pcmtx,HH,pcoptions);
end

[v1,posdef,pcgit] = pcgr(DM,DG,grad,kmax,tol,mtxmpy,data,H,R,permR);

% FORM A 2-D SUBSPACE
v1 = v1/norm(v1);
Z(:,1) = v1;
if (posdef < 1)
   v2 = D*sign(grad);
   v2 = v2/norm(v2);
   v2 = v2 - v1*(v1'*v2);
   nrmv2 = norm(v2);
   if nrmv2 > tol2
      v2 = v2/nrmv2; Z(:,2) = v2;
   end
else
   v1 = v1/norm(v1); 
   v2 = grad/norm(grad);
   v2 = v2 - v1*(v1'*v2); 
   nrmv2 = norm(v2);
   if nrmv2 > tol2
      v2 = v2/nrmv2; Z(:,2) = v2;
   end
end

W = DM*Z; 
WW = feval(mtxmpy,W,H,data); 
W = DM*WW;
MM = Z'*W + Z'*DG*Z;
rhs=full(Z'*grad);

% SOLVE TRUST REGION OVER 2-D.
[ss,qpval,po,fcnt,lambda] = trust(rhs,MM,delta);
ss = Z*ss;
norms = norm(ss);
s = abs(diag(D)).*ss;




