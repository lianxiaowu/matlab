function[x,val,output,ex,lambda]=sqpmin(c,H,xstart,A,b,lb,ub,verb,options,computeLambda)
%SQPMIN	Solve quadratic problems with box constraints or linear equalities
%
% Locate local soln to
%
%        min { q(x) = .5x'Hx + c'x : l <= x <= u}. 
%
%                            or
%
%       min { q(x) = .5x'Hx + c'x : Ax = b},
%
%
% where H is sparse symmetric mtx. (may be virtual),
%
% x = sqpmin(c,H,xstart,options) return the minimizer of the
% quadratic function q(x) subject to any bounds indicated in
% the named parameter list options. xstart is the starting point.
%
% x = sqpmin(c,H,xstart,options,A,b) solves the linearly
% constrained problem, min { q(x) = .5x'Hx + c'x : Ax = b}.
%
% [x,val] =  sqpmin(c,H,xstart,A,b,ub,lb) returns the value of 
% the quadratic objective function at the solution.
%
% [x,val,gopt] = sqpmin(c,H,xstart,options, ...) returns a measure
% of first-order optimality.
%
% [x,val,gopt,it] = sqpmin(c,H,xstart,options,...) returns
% number of iterations used.
%
% [x,val,gopt,it,npcg] =  sqpmin(c,H,xstart,options,...) returns
% total number of conjugate gradient iterations used.
%
% [x,val,gopt,it,npcg,ex] = sqpmin(c,H,xstart,options,...) returns
% termination code.

%   Copyright (c) 1990-98 by The MathWorks, Inc.
%   $Revision: 1.13 $  $Date: 1998/09/09 19:48:01 $

if nargin < 2, error('sqpmin requires at least 2 arguments'), end
if nargin <=2, xstart = []; end
n = length(c); 
H = sparse(H);

if isempty(lb), lb = -inf*ones(n,1); end
if isempty(ub),ub = inf*ones(n,1); end
arg = (ub >= 1e10); arg2 = (lb <= -1e10);
ub(arg) = inf;
lb(arg2) = -inf;
if any(ub == lb) 
   errmsg=sprintf('%s\n%s',...
      'Equal upper and lower bounds not permitted in this large-scale method.',...
      'Use equality constraints and the medium-scale method instead.');
   error(errmsg)
elseif min(ub-lb) <= 0
   error('Inconsistent bounds.')
end
if isempty(xstart), xstart = startx(ub,lb); end
if min(min(ub-xstart),min(xstart-lb)) < 0, xstart = startx(ub,lb); end

% get options out
typx = optimget(options,'typicalx') ;
% In case the defaults were gathered from calling: optimset('quadprog'):
numberOfVariables = n;
if ischar(typx)
   typx = eval(typx);
end

% Later ShowStatusWindow, Preconditioner and HessMult will be user-settable
pcmtx = optimget(options,'Preconditioner','hprecon') ;
mtx = optimget(options,'HessMult','hmult') ;
showstat = optimget(options,'showstatus','off');
kmax = optimget(options,'MaxPCGIter', max(1,floor(n/2))) ;
if ischar(kmax)
   kmax = eval(kmax);
end

pcf = optimget(options,'PrecondBandWidth') ;
tolx = optimget(options,'tolx') ;
tolfun = optimget(options,'tolfun');
itb = optimget(options,'maxiter') ;
switch showstat
case 'iter'
   showstat = 2;
case {'none','off'}
   showstat = 0;
case 'final'
   showstat = 1;
case 'iterplus'
   showstat = 3;
otherwise
   showstat = 0;
end

 
if n == 0, error('n must be positive'), end

if ~showstat, delete(findobj('type','figure',...
      'Name','Algorithm Performance Statistics')) ;
end ;


%   INITIALIZATIONS
lambda.lower = [];
lambda.upper = [];
lambda.eqlin = [];  
lambda.ineqlin = [];  % This won't change because no inequalities.
val = []; gopt=[];
output = [];
it = 1; cvec = c; nbnds = 1;
if nargin < 5, A = []; end
if isempty(A)
   fdata = [];
   
   % Box-constrained problem
   
   [x,val,gopt,it,npcg,ex,lambda]=sqpbox(c,H,lb,ub,xstart,typx,verb,pcmtx,pcf,...
     mtx,fdata,tolx,tolfun,itb,showstat,computeLambda,kmax);
  lambda.ineqlin = []; lambda.eqlin = [];
  output.firstorderopt = gopt;
  output.iterations = it; 
  output.cgiterations = npcg;
  output.algorithm = 'large-scale: reflective trust-region';
else
   if ((max(lb) > -inf) | (min(ub) < inf))
      error('sqpmin doesn''t handle both box constraints and Ax = b');
   else
      
      %          Equality constrained problem
      [mA,nA] = size(A);
      if nargin < 6, b = zeros(mA,1); end
      if isempty(b), b = zeros(mA,1); end
      fdata = [];
      % Note we pass options in so some values are different for PPCGR than for SQPBOX
      [x,po,npcg,pgnrm,ex,lambda]=ppcgr(c,H,A,b,options,verb,fdata,computeLambda);
      if ex == -2
         % ppcgr aborted
         return
      end
      
      it = 1; 
      w = feval(mtx,x,H,fdata); 
      g = c + w;
      gopt = pgnrm;
      val = x'*(c + .5*w);
      output.firstorderopt = gopt;
      output.iterations = it; 
      output.cgiterations = npcg;
      output.algorithm = 'large-scale: projective preconditioned conjugate gradients';
   end
end

% Temporary until we make the other flag values consistent
if ex == 4
   ex = 0;
elseif ex > 0
   ex =1;
end




