 function[x,FVAL,LAMBDA,EXITFLAG,OUTPUT,GRAD,HESSIAN]=sfminbx(funfcn,x,l,u,verb,options,...
   computeLambda,initialf,initialGRAD,initialHESS,Hstr,varargin)
%SFMINBX Nonlinear minimization with box constraints.
%
% [x,val,gopt,it,npcg,ex]=sfminbx(fname,xstart,fdata,l,u,verb,...
%	             pcmtx,pcflags,mtxmpy,tol,maxiter,showstat,Hstr)
% Locate a local minimizer to 
%
%               min { f(x) :  l <= x <= u}.
%
%	where f(x) maps n-vectors to scalars.
% The driver function is SFMIN

%   Copyright (c) 1990-98 by The MathWorks, Inc.
%   $Revision: 1.20 $  $Date: 1998/09/15 21:13:36 $

%   Initialization
xcurr = x(:);  % x has "the" shape; xcurr is a vector
n = length(xcurr); 
it= 1; totm = 0; totls = 0;

header = sprintf(['\n                                Norm of      First-order \n',...
      ' Iteration        f(x)          step          optimality   CG-iterations']);
formatstr = ' %5.0f      %13.6g  %13.6g   %12.3g     %7.0f';
if n == 0, 
   error('n must be positive'), 
end
fdata = [];
if isempty(l), 
   l = -inf*ones(n,1); 
end, 
if isempty(u), 
   u = inf*ones(n,1);
end
arg = (u >= 1e10); arg2 = (l <= -1e10);
u(arg) = inf*ones(length(arg(arg>0)),1);
l(arg2) = -inf*ones(length(arg2(arg2>0)),1);
if any(u == l) 
   errmsg=sprintf('%s\n%s',...
      'Equal upper and lower bounds not permitted in this large-scale method.',...
      'Use equality constraints and the medium-scale method instead.');
   error(errmsg)
elseif min(u-l) <= 0
   error('Inconsistent bounds.')
end
if min(min(u-xcurr),min(xcurr-l)) < 0, xcurr = startx(u,l); end

% get options out
typx = optimget(options,'typicalx') ;
% In case the defaults were gathered from calling: optimset('quadprog'):
numberOfVariables = n;
if ischar(typx)
   typx = eval(typx);
end

% Will be user-settable later:
showstat = optimget(options,'showstatus','off');
pcmtx = optimget(options,'Preconditioner','hprecon') ;
mtxmpy = optimget(options,'HessMult','hmult') ;

switch showstat
case 'iter'
   showstat = 2;
case {'none','off'}
   showstat = 0;
case 'final'
   showstat = 1;
case {'iterplus','iterplusbounds'}  % if no finite bounds, this is the same as 'iter'
   showstat = 3;
otherwise
   showstat = 1;
end
active_tol = optimget(options,'ActiveConstrTol',sqrt(eps));
pcflags = optimget(options,'PrecondBandWidth') ;
tol2 = optimget(options,'tolx') ;
tol1 = optimget(options,'tolfun') ;
tol = tol1;
maxiter = optimget(options,'maxiter') ;
maxfunevals = optimget(options,'maxFunEvals') ;
pcgtol = optimget(options,'TolPCG', 0.1) ;  % pcgtol = .1;
kmax = optimget(options,'MaxPCGIter', max(1,floor(n/2))) ;
if ischar(kmax)
   kmax = eval(kmax);
end
if ischar(maxfunevals)
   maxfunevals = eval(maxfunevals);
end
maxcount = min(maxiter, maxfunevals); % numfunevals = iterations, so just take minimum

dnewt = []; gopt = [];
ex = 0; posdef = 1; npcg = 0; 

%tol1 = tol; tol2 = sqrt(tol1)/10; 
if strcmp(optimget(options,'DerivativeCheck'),'on')
   warnstr = sprintf('%s\n%s\n', ...
      'Trust region algorithm does not currently check user-supplied gradients,', ...
      '  ignoring OPTIONS.DerivativeCheck.');
   warning(warnstr);
end

vpos(1,1) = 1; vpcg(1,1) = 0; nbnds = 1;
pcgit = 0; delta = 10;nrmsx = 1; ratio = 0; degen = inf; 
if (all(u == inf) & all(l == -inf)) nbnds = 0; end
DS = speye(n);   v = zeros(n,1); dv = ones(n,1); del = 10*eps;
oval = inf;  g = zeros(n,1); newgrad = g; Z = []; 

% Make x conform to the user's input x
x(:) = xcurr;
%   Evaluate f,g,  and H
if ~isempty(Hstr)  % use sparse finite differencing
   %[val,g] = feval(fname,x,fdata);
   switch funfcn{1}
   case 'fun'
      error('should not reach this')
   case 'fungrad'
      %[val,g(:)] = feval(funfcn{3},x,varargin{:});
      val = initialf; g(:) = initialGRAD;
   case 'fun_then_grad'
      % val = feval(funfcn{3},x,varargin{:}); 
      % g(:) = feval(funfcn{4},x,varargin{:});
      val = initialf; g(:) = initialGRAD;
   otherwise
      if isequal(funfcn{2},'fmincon')
         error('Undefined calltype in FMINCON');
      else
         error('Undefined calltype in FMINUNC');
      end
   end
   
   %      Determine coloring/grouping for sparse finite-differencing
   p = colmmd(Hstr)'; p = (n+1)*ones(n,1)-p; group = color(Hstr,p);
   % pass in the user shaped x
   H = sfd(x,g,Hstr,group,[],funfcn,varargin{:});
   %
else % user-supplied computation of H or dnewt
   % [val,g,H] = feval(fname,x,fdata);
   switch funfcn{1}
   case 'fungradhess'
     % [val,g(:),H] = feval(funfcn{3},x,varargin{:});
      val = initialf; g(:) = initialGRAD; H = initialHESS;
   case 'fun_then_grad_then_hess'
      % val = feval(funfcn{3},x,varargin{:}); 
      % g(:) = feval(funfcn{4},x,varargin{:});
      % H = feval(funfcn{5},x,varargin{:});  
      val = initialf; g(:) = initialGRAD; H = initialHESS;
   otherwise
      if isequal(funfcn{2},'fmincon')
         error('Undefined calltype in FMINCON');
      else
         error('Undefined calltype in FMINUNC');
      end
   end
end

delbnd = max(100*norm(xcurr),1);
[nn,pp] = size(g);

%   Extract the Newton direction?
if pp == 2, dnewt = g(1:n,2); end
if showstat > 1
   figtr=display1('init',maxiter,tol,showstat,nbnds,xcurr,g(:,1),l,u); 
end
if verb > 1
   disp(header)
end

%   MAIN LOOP: GENERATE FEAS. SEQ.  xcurr(it) S.T. f(x(it)) IS DECREASING.
while ~ex
   if ~isfinite(val) | any(~isfinite(g))
      errmsg= sprintf('%s%s%s',funfcn{2},' cannot continue: ',...
         'user function is returning Inf or NaN values.');
      error(errmsg)
   end
   
   %     Stop (interactive)?
   figtr = findobj('type','figure','Name','Progress Information') ;
   if ~isempty(figtr)
      lsotframe = findobj(figtr,'type','uicontrol',...
         'Userdata','LSOT frame') ;
      if get(lsotframe,'Value'), 
         ex = 10 % New exiting condition 
         if verb > 0
            display('Exiting per request.')
         end
      end 
   end 
   
   %     Update
   [v,dv] = definev(g(:,1),xcurr,l,u); 
   gopt = v.*g(:,1); gnrm = norm(gopt,inf);
   vgnrm(it,1)=gnrm;
   r = abs(min(u-xcurr,xcurr-l)); degen = min(r + abs(g(:,1)));
   vdeg(it,1) = min(degen,1); bndfeas = min(min(xcurr-l,u-xcurr));
   if ((u == inf*ones(n,1)) & (l == -inf*ones(n,1))) degen = -1; end
   
   % Display
   if showstat > 1
      display1('progress',it,gnrm,val,pcgit,npcg,degen,...
         bndfeas,showstat,nbnds,xcurr,g(:,1),l,u,figtr,posdef); 
   end
   if verb > 1
      currOutput = sprintf(formatstr,it,val,nrmsx,gnrm,pcgit);
      disp(currOutput);
   end
   
   %     TEST FOR CONVERGENCE
   diff = abs(oval-val); 
   oval = val; vflops(it,1) = flops; totm = flops;
   if (nrmsx < .9*delta)&(ratio > .25)&(diff < tol1*(1+abs(oval)))
      ex = 1;
      if verb > 0
         disp('Optimization terminated successfully:')
         disp(' Relative function value changing by less than OPTIONS.TolFun');
      end
      
   elseif (it > 1) & (nrmsx < tol2) 
      ex = 2;
      if verb > 0
         disp('Optimization terminated successfully:')
         disp(' Norm of the current step is less than OPTIONS.TolX');
      end
      
   elseif ((gnrm < tol1) & (posdef ==1) )
      ex = 3;
      if verb > 0
         disp('Optimization terminated successfully:')
         disp(' First-order optimality less than OPTIONS.TolFun, and no negative/zero curvature detected');
      end
   end
   
   %     Step computation
   if ~ex
      
      %       Determine trust region correction
      dd = abs(v); D = sparse(1:n,1:n,full(sqrt(dd))); 
      sx = zeros(n,1); theta = max(.95,1-gnrm);  
      oposdef = posdef;
      [sx,snod,qp,posdef,pcgit,Z] = trdog(xcurr, g(:,1),H,fdata,D,delta,dv,...
         mtxmpy,pcmtx,pcflags,pcgtol,kmax,theta,l,u,Z,dnewt,'hessprecon');
      if isempty(posdef), posdef = oposdef; end
      nrmsx=norm(snod); npcg=npcg + pcgit;
      newx=xcurr + sx; vpcg(it+1,1)=pcgit;
      
      %       Perturb?
      [pert,newx] = perturb(newx,l,u);
      vpos(it+1,1) = posdef; 
      
      % Make newx conform to user's input x
      x(:) = newx;
      % Evaluate f, g, and H
      if ~isempty(Hstr) % use sparse finite differencing
         %[newval,newgrad] = feval(fname,x,fdata);
         switch funfcn{1}
         case 'fun'
            error('should not reach this')
         case 'fungrad'
            [newval,newgrad(:)] = feval(funfcn{3},x,varargin{:});
         case 'fun_then_grad'
            newval = feval(funfcn{3},x,varargin{:}); 
            newgrad(:) = feval(funfcn{4},x,varargin{:});
         otherwise
            error('Undefined calltype in FMINUNC');
         end
         newH = sfd(x,newgrad,Hstr,group,[],funfcn,varargin{:});
         
      else % user-supplied computation of H or dnewt
         %[newval,newgrad,newH] = feval(fname,x,fdata);
         switch funfcn{1}
         case 'fungradhess'
            [newval,newgrad(:),newH] = feval(funfcn{3},x,varargin{:});
         case 'fun_then_grad_then_hess'
            newval = feval(funfcn{3},x,varargin{:}); 
            newgrad(:) = feval(funfcn{4},x,varargin{:});
            newH = feval(funfcn{5},x,varargin{:});
         otherwise
            error('Undefined calltype in FMINUNC');
         end
         
      end
      [nn,pp] = size(newgrad);
      aug = .5*snod'*((dv.*abs(newgrad(:,1))).*snod);
      ratio = (newval + aug -val)/qp; vratio(it,1) = ratio;
      
      if (ratio >= .75) & (nrmsx >= .9*delta)
         delta = min(delbnd,2*delta);
      elseif ratio <= .25
         delta = min(nrmsx/4,delta/4);
      end
      if newval == inf
         delta = min(nrmsx/20,delta/20);
      end
      
      %       Update
      if newval < val
         xold = xcurr; xcurr=newx; val = newval; g= newgrad; H = newH;
         Z = [];
         
         %          Extract the Newton direction?
         if pp == 2, dnewt = newgrad(1:n,2); end
      end
      it = it+1; vval(it,1) = val;
   end
   if it > maxcount, 
      ex=4; 
      it = it-1; 
      if verb > 0
         if it > maxiter
            disp('Maximum number of iterations exceeded;')
            disp('   increase options.MaxIter')
         elseif it > maxfunevals
            disp('Maximum number of function evaluations exceeded;')
            disp('   increase options.MaxFunEvals')
         end
      end
   end              
end % while

if showstat > 1
   display1('final',figtr); 
end
if showstat 
   xplot(it,vval,vgnrm,vflops,vpos,vdeg,vpcg);
end
HESSIAN = H;
GRAD = g;
FVAL = val;
LAMBDA = [];
if ex==4
   EXITFLAG = 0;
elseif ex==10
   EXITFLAG = -1;
else
   EXITFLAG = 1;
end
OUTPUT.iterations = it;
OUTPUT.funcCount = it;
OUTPUT.cgiterations = npcg;
OUTPUT.firstorderopt = gnrm;
OUTPUT.algorithm = 'large-scale: trust-region reflective Newton'; 
x(:) = xcurr;
if computeLambda
   g = full(g);
   
   LAMBDA.lower = zeros(length(l),1);
   LAMBDA.upper = zeros(length(u),1);
   argl = logical(abs(xcurr-l) < active_tol);
   argu = logical(abs(xcurr-u) < active_tol);
   
   LAMBDA.lower(argl) = (g(argl));
   LAMBDA.upper(argu) = -(g(argu));
   LAMBDA.ineqlin = []; LAMBDA.eqlin = []; LAMBDA.ineqnonlin=[]; LAMBDA.eqnonlin=[];
else
   LAMBDA = [];   
end




