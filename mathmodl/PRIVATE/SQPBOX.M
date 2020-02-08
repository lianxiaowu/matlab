function[x,val,csnrm,it,npcg,ex,LAMBDA]=sqpbox(c,H,l,u,xstart,typx,verb,...
   pcmtx,pcflags,mtxmpy,data,tolx,tolfun,itb,showstat,computeLambda,kmax)
%SQPBOX Minimize box-constrained quadratic fcn
%
%   [x,val,csnrm,it,npcg,ex]=sqpbox(c,H,l,u,xstart,typx,verb,...
%                   pcmtx,pcflags,mtxmpy,data,tol,itb,showstat)
%
%   Locate a (local) soln
%   to the box-constrained QP:
%
%        min { q(x) = .5x'Hx + c'x :  l <= x <= u}. 
%
%   where H is sparse symmetric, c is a col vector,
%   l,u are vectors of lower and upper bounds respectively.
% Driver function is SQPMIN

%   Copyright (c) 1990-98 by The MathWorks, Inc.
%   $Revision: 1.8 $  $Date: 1998/08/17 19:42:12 $

%   INITIALIZATIONS
if nargin <= 1
   error('sqpbox requires at least 2 arguments');
end
n = length(c); 
it = 0; 
cvec = c; 
nbnds = 1;
header = sprintf(['\n                                Norm of      First-order \n',...
      ' Iteration        f(x)          step          optimality   CG-iterations']);
formatstr = ' %5.0f      %13.6g  %13.6g   %12.3g     %7.0f';

if n == 0
   error('n must be positive')
end
if nargin <= 2
   l = -inf*ones(n,1);
end
if nargin <= 3, 
   u = inf*ones(n,1); 
end
if isempty(l),
   l = -inf*ones(n,1); 
end 
if isempty(u),
   u = inf*ones(n,1);
end
arg = (u >= 1e10); 
arg2 = (l <= -1e10);
u(arg) = inf;
l(arg2) = -inf;
if min(u-l) <= 0, 
   error('inconsistent bounds'),
end
lvec = l; uvec = u; 
if nargin <= 4, 
   xstart = startx(u,l);
end
if min(min(u-xstart),min(xstart-l)) < 0, 
   xstart = startx(u,l);
end
if nargin <=5, 
   typx = ones(n,1);
end
if isempty(typx), 
   typx = ones(n,1);
end
if nargin <=6, 
   verb = 0; 
end
if isempty(verb), 
   verb = 0; 
end
if nargin <= 7, 
   pcmtx = 'hprecon', 
end
if isempty(pcmtx), 
   pcmtx = 'hprecon', 
end
if nargin <= 8, 
   pcflags = [0;0]; 
end
if isempty(pcflags), 
   pcflags = [0;0]; 
end
if nargin <= 9, 
   mtxmpy = 'hmult',
end
if isempty(mtxmpy), 
   mtxmpy = 'hmult', 
end
if nargin <= 10, 
   data = [];
end
if nargin <= 11, 
   tolx = 100*eps; 
end
if nargin <= 12
   tolfun = 100*eps;
end
if isempty(tolx), 
   tolx = (10^2)*eps; 
end
if isempty(tolfun)
   tolfun = 100*eps;
end
if nargin <=13, 
   itb = 200; 
end
if isempty(itb), 
   itb = 200 ; 
end  
if nargin <=14, 
   showstat = 0; 
end
if isempty(showstat), 
   showstat = 0 ; 
end

pcgit = 0;
pcgtol = .1; 
tolx2 = sqrt(tolx); 
tolfun2 = sqrt(tolfun);
vpcg(1,1) = 0; 
vpos(1,1) = 1;
[xstart,l,u,ds,DS,c] = shiftsc(xstart,l,u,typx,mtxmpy,data,cvec,H);
dellow = 1.; 
delup = 10^3; 
npcg = 0; 
digits = inf; 
ex = 0; 
v = zeros(n,1);
dv = ones(n,1);
del = 10*eps;
posdef = 1;
x = xstart; 
y = x;
sigma = ones(n,1);
g = zeros(n,1);
oval = inf; 
[val,g] = fquad(x,c,H,mtxmpy,data,DS); 

if ((u == inf*ones(n,1)) & (l == -inf*ones(n,1))) 
   nbnds = 0; 
end
if showstat > 1, 
   figtr = display1('init',itb,tolfun,showstat,nbnds,x,g,l,u); 
end
%
%   MAIN LOOP: GENERATE FEAS. SEQ.  x(it) S.T. q(x(it)) IS DECREASING.
while ~ex
   it = it + 1;
   vval(it,1) = val;
   %
   %     terminate?
   figtr = findobj('type','figure','Name','Progress Information') ;
   if ~isempty(figtr)
      lsotframe = findobj(figtr,'type','uicontrol',...
         'Userdata','LSOT frame') ;
      if get(lsotframe,'Value'), 
         ex = 10; % New exiting condition  
      end ;
   end ;
   %
   %     Update and display
   [v,dv] = definev(g,x,l,u); 
   csnrm = norm(v.*g,inf); 
   vcsnrm(it,1) = csnrm;
   r = abs(min(u-x,x-l));
   degen = min(r + abs(g));
   vdeg(it,1) = min(degen,1);
   if ((u == inf*ones(n,1)) & (l == -inf*ones(n,1))) 
      degen = -1; 
   end
   bndfeas = min(min(x-l,u-x));
   if showstat > 1
      display1('progress',it,csnrm,val,pcgit,npcg,degen,...
         bndfeas,showstat,nbnds,x,g,l,u,figtr);
   end
   %
   %     TEST FOR CONVERGENCE
   diff = abs(oval-val); 
   vdiff(it,1) = diff/(1+abs(oval));
   if it > 1, 
      digits = (prev_diff)/max(diff,eps); 
   end
   prev_diff = diff; 
   oval = val; 
   vflops(it,1) = flops;
   if diff < tolfun*(1+abs(oval)),
      ex = 1;
      if verb > 0
         disp('Optimization terminated successfully:')
         disp(' Relative function value changing by less than OPTIONS.TolFun');
      end
      
   elseif ((diff < tolfun2*(1+abs(oval))) & (digits < 3.5)) & posdef, 
      ex = 2;
      if verb > 0
         disp('Optimization terminated successfully:')
         disp(' Relative function value changing by less than sqrt(OPTIONS.TolFun), ');
         disp(' no negative curvature detected in Hessian this iteration, and ');
         disp(' the rate of progress (change in f(x)) is slow');
      end
   elseif ((csnrm < tolfun) & posdef), 
      ex = 3;
      if verb > 0
         disp('Optimization terminated successfully:')
         disp(' No negative curvature detected in Hessian this iteration, and ');
         disp(' first order optimality measure < OPTIONS.TolFun');
      end
   end
   %
   if ~ex 
      %       DETERMINE THE SEARCH DIRECTION
      dd = abs(v); 
      D = sparse(1:n,1:n,full(sqrt(dd).*sigma));
      grad = D*g; 
      normg = norm(grad); 
      
      delta = max(dellow,norm(v)); 
      delta = min(delta,delup); 
      vdelta(it,1) = delta;
      [s,posdef,pcgit] = drqpbox(D,DS,grad,delta,g,dv,mtxmpy,data,...
         pcmtx,pcflags,pcgtol,H,0,kmax);
      npcg = npcg + pcgit; 
      vpos(it+1,1) = posdef; 
      vpcg(it+1,1) = pcgit;
      %
      %       DO A REFLECTIVE (BISECTION) LINE SEARCH. UPDATE x,y,sigma.
      strg= s'*(sigma.*g); 
      ox = x; 
      osig = sigma; 
      ostrg = strg;
      if strg >= 0, 
         ex = -1; %ex = 5; 
         if verb > 0
            disp('Optimization terminated:')
            disp(' Loss of feasibility with respect to the constraints detected.');
         end
         
      else,
         [x,sigma,alpha] = biqpbox(s,c,ostrg,ox,y,osig,l,u,oval,posdef,...
            normg,DS,mtxmpy,data,H);
         if alpha == 0, 
            ex = 6; 
            if verb > 0
               disp('Optimization terminated:')
               disp(' Current direction not descent direction; ')
               disp(' the problem may be ill-conditioned.');
            end
            
         end
         y = y + alpha*s; 
         %
         %          PERTURB x AND y ?
         [pert,x,y] = perturb(x,l,u,del,y,sigma);
         %
         %          EVALUATE NEW FUNCTION VALUE, GRADIENT. 
         [val,g] = fquad(x,c,H,mtxmpy,data,DS); 
      end
      if it >= itb, 
         ex=4; 
         if verb > 0
            disp('Maximum number of iterations exceeded;')
            disp('   increase options.MaxIter')
         end              
      end
   end
end
%
%   RESCALE, UNSHIFT, AND EXIT.
x = unshsca(x,lvec,uvec,DS);
[val,g] = fquad(x,cvec,H,mtxmpy,data); 
if showstat > 1,
   display1('final',figtr);
end
if showstat>0, 
   xplot(it,vval,vcsnrm,vflops,vpos,vdeg,vpcg); 
end

if computeLambda
 LAMBDA.lower = zeros(length(lvec),1);
 LAMBDA.upper = zeros(length(uvec),1);
 active_tol = sqrt(eps);
 argl = logical(abs(x-lvec) < active_tol);
 argu = logical(abs(x-uvec) < active_tol);

 g = full(g);
 LAMBDA.lower(argl) = abs(g(argl));
 LAMBDA.upper(argu) = -abs(g(argu));
else
  LAMBDA=[];
end
