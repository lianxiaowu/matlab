function[x,FVAL,LAMBDA,EXITFLAG,OUTPUT,GRAD,HESSIAN]= sfminle(funfcn,x,A,b,verb,options,...
   computeLambda,initialf,initialGRAD,initialHESS,Hstr,varargin)
%SFMINLE Nonlinear minimization with linear equalities.
%
% [x,val,g,it,npcg,ex]=sfminle(fname,xstart,A,b,fdata,verb,...
%	         pcmtx,pcflags,mtxmpy,tol,itb,showstat,Hstr,options)
% Locate a local minimizer to 
%
%               min { f(x) :  Ax = b}.
%
%	where f(x) maps n-vectors to scalars.
%
% Driver function is SFMIN

%   Copyright (c) 1990-98 by The MathWorks, Inc.
%   $Revision: 1.14 $  $Date: 1998/09/15 21:13:37 $

%   Initialization
xcurr = x(:); % x has "the" shape; xcurr is a vector

n = length(xcurr); it= 1; totm = 0; totls = 0;
header = sprintf(['\n                                Norm of      First-order \n',...
      ' Iteration        f(x)          step          optimality   CG-iterations']);
formatstr = ' %5.0f      %13.6g  %13.6g   %12.3g     %7.0f';

if n == 0, 
   error('n must be positive'), 
end
fdata = [];
[mm,nn] = size(A);
if n ~= nn
   error('Column dimension of Aeq is inconsistent with length of x'), end
m = length(b);
if m ~= mm
   error('Row dimension of Aeq is inconsistent with length of beq'), end

% get options out
typx = optimget(options,'typicalx') ;
% In case the defaults were gathered from calling: optimset('quadprog'):
numberOfVariables = n;
if ischar(typx)
   typx = eval(typx);
end
% Will be user-settable later:
showstat = optimget(options,'showstatus','off');
pcmtx = optimget(options,'Preconditioner','preaug') ;
if isequal(pcmtx,'hprecon')
   error('Preconditioner for equality constrained problem should be preaug')
end
mtxmpy = optimget(options,'HessMult','hmult') ;

switch showstat
case 'iter'
   showstat = 2;
case {'none','off'}
   showstat = 0;
case 'final'
   showstat = 1;
case 'iterplusbounds'  % if no finite bounds, this is the same as 'iter'
   showstat = 3;
otherwise
   showstat = 1;
end

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
% tol = (10^6)*eps; 
%tol1 = tol; tol2 = sqrt(tol1)/10; 

if strcmp(optimget(options,'DerivativeCheck'),'on')
   warnstr = sprintf('%s\n%s\n', ...
      'Trust region algorithm does not currently check user-supplied gradients,', ...
      '  ignoring OPTIONS.DerivativeCheck.');
   warning(warnstr);
end

dnewt = []; gopt = []; nrows = 0;
snod = [];
ex = 0; posdef = 1; npcg = 0; vpos(1,1) = 1; vpcg(1,1) = 0; 
pcgit = 0; delta = 100;nrmsx = 1; ratio = 0; degen = inf; 
pcgtol = 1e-2;
DS = speye(n);   v = ones(n,1); dv = ones(n,1); del = 10*eps;
oval = inf;  gradf = zeros(n,1); newgrad = gradf; Z = []; 

%   Remove (numerical) linear dependencies from A
AA = A; bb = b;
[A,b] = dep(AA,[],bb);
[m,n1] = size(A); [mm,n2] = size(AA);
if verb > 1 & m ~= mm
   disp('linear dependencies'), mm-m, disp(' rows of A removed'); end

%   Get feasible: nearest feas. pt. to xstart
xcurr = feasibl(A,b,xcurr);


% Make x conform to the user's input x
x(:) = xcurr;

%   Evaluate f,g,  and H
if ~isempty(Hstr) % use sparse finite differencing
   % [val,gradf] = feval(fname,x,fdata);
         
      switch funfcn{1}
      case 'fun'
         error('should not reach this')
      case 'fungrad'
         val = initialf; gradf(:) = initialGRAD;
        % [val,gradf(:)] = feval(funfcn{3},x,varargin{:});
     case 'fun_then_grad'
        val = initialf; gradf(:) = initialGRAD;
        %   val = feval(funfcn{3},x,varargin{:}); 
        %   gradf(:) = feval(funfcn{4},x,varargin{:});
     otherwise
        error('Undefined calltype in FMINCON');
      end
      
   %      Determine coloring/grouping for sparse finite-differencing
   p = colmmd(Hstr)'; p = (n+1)*ones(n,1)-p; group = color(Hstr,p);
   H = sfd(x,gradf,Hstr,group,[],funfcn,varargin{:});

else % user-supplied computation of H or dnewt
   %[val,gradf,H] = feval(fname,x,fdata);
   switch funfcn{1}
   case 'fungradhess'
      val = initialf; gradf(:) = initialGRAD; H = initialHESS;
      % [val,gradf(:),H] = feval(funfcn{3},x,varargin{:});
   case 'fun_then_grad_then_hess'
       val = initialf; gradf(:) = initialGRAD; H = initialHESS;
      % val = feval(funfcn{3},x,varargin{:}); 
      % gradf(:) = feval(funfcn{4},x,varargin{:});
      % H = feval(funfcn{5},x,varargin{:});
     
   otherwise
      error('Undefined calltype in FMINCON');
   end
end
[nn,pnewt] = size(gradf);

%   Extract the Newton direction?
if pnewt == 2, dnewt = gradf(1:n,2); end
PT = findp(A);
[g,LZ,UZ,pcolZ,PZ] = project(A,-gradf(1:n,1),PT);
gnrm = norm(g,inf);

if showstat > 1, 
   figtr=display1('init',maxiter,tol,showstat,0,xcurr,g(:,1),[],[]); 
end
if verb > 1
   disp(header)
end

%   MAIN LOOP: GENERATE FEAS. SEQ.  xcurr(it) S.T. f(xcurr(it)) IS DECREASING.
while ~ex
   if ~isfinite(val) | any(~isfinite(gradf))
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
         EXITFLAG = -1;
         if verb > 0
            display('Exiting per request.')
         end

      end 
   end 

   %     Update and display
   vgnrm(it,1)=gnrm;
   if showstat > 1
      display1('progress',it,gnrm,val,pcgit,npcg,degen,...
         [],showstat,0,xcurr,g(:,1),[],[],figtr,posdef); 
   end
   if verb > 1
      currOutput = sprintf(formatstr,it,val,nrmsx,gnrm,pcgit);
      disp(currOutput);
   end

   %     TEST FOR CONVERGENCE
   diff = abs(oval-val); 
   oval = val; vflops(it,1) = flops; totm = flops;
   if (nrmsx < .9*delta)&(ratio > .25)&(diff < tol1*(1+abs(oval)))
      ex = 1; EXITFLAG = 1;
      if verb > 0
         disp('Optimization terminated successfully:')
         disp(' Relative function value changing by less than OPTIONS.TolFun');
      end

elseif (it > 1) & (nrmsx < tol2), 
   ex = 2; EXITFLAG = 1;
   if verb > 0
         disp('Optimization terminated successfully:')
         disp(' Norm of the current step is less than OPTIONS.TolX');
      end

elseif ((gnrm < tol1) & (posdef ==1) ),
   ex = 3; EXITFLAG = 1;
   if verb > 0
         disp('Optimization terminated successfully:')
         disp(' First-order optimality less than OPTIONS.TolFun, and no negative/zero curvature detected');
      end
   end

   %     Step computation
   if ~ex

      %       Determine trust region correction
      dd = abs(v); D = sparse(1:n,1:n,full(sqrt(dd))); 
      grad = D*g(:,1);
      sx = zeros(n,1); theta = max(.95,1-gnrm);  
      oposdef = posdef;
      [sx,snod,qp,posdef,pcgit,Z] = trdg(xcurr,gradf(:,1),H,fdata,...
         delta,g,mtxmpy,pcmtx,pcflags,...
         pcgtol,kmax,A,zeros(m,1),Z,dnewt,options,...
         PT,LZ,UZ,pcolZ,PZ);

if isempty(posdef), 
   posdef = oposdef; 
end
nrmsx=norm(snod); 
npcg=npcg + pcgit;
newx=xcurr + sx; 
vpcg(it+1,1)=pcgit;
      vpos(it+1,1) = posdef;
      
      % Make newx conform to user's input x
      x(:) = newx;
      %       Evaluate f,g,  and H
      if ~isempty(Hstr) % use sparse finite differencing
         % [newval,newgrad] = feval(fname,newx,fdata);
         switch funfcn{1}
         case 'fun'
            error('should not reach this')
         case 'fungrad'
            [newval,newgrad(:)] = feval(funfcn{3},x,varargin{:});
            %OPTIONS(11)=OPTIONS(11)+1;
         case 'fun_then_grad'
            newval = feval(funfcn{3},x,varargin{:}); 
            newgrad(:) = feval(funfcn{4},x,varargin{:});
            % OPTIONS(11)=OPTIONS(11)+1;      
         otherwise
            error('Undefined calltype in FMINUNC');
         end

         newH = sfd(x,newgrad,Hstr,group,[],funfcn,varargin{:});

      else % user-supplied computation of H or dnewt
         %[newval,newgrad,newH] = feval(fname,newx,fdata);
         switch funfcn{1}
         case 'fungradhess'
            [newval,newgrad(:),newH] = feval(funfcn{3},x,varargin{:});
            %OPTIONS(11)=OPTIONS(11)+1;
         case 'fun_then_grad_then_hess'
            newval = feval(funfcn{3},x,varargin{:}); 
            newgrad(:) = feval(funfcn{4},x,varargin{:});
            newH = feval(funfcn{5},x,varargin{:});
            % OPTIONS(11)=OPTIONS(11)+1;
            
         otherwise
            error('Undefined calltype in FMINUNC');
         end

      end
      [nn,pnewt] = size(newgrad);
      if pnewt == 2, 
         dnewt = newgrad(1:n,2); 
      end
      aug = .5*snod'*((dv.*abs(newgrad(:,1))).*snod);
      ratio = (newval + aug -val)/qp; 
      vratio(it,1) = ratio;
      if (ratio >= .75) & (nrmsx >= .9*delta),
         delta = 2*delta;
      elseif ratio <= .25, 
         delta = min(nrmsx/4,delta/4);
      end
      if newval == inf; 
         delta = min(nrmsx/20,delta/20);
      end

      %       Update
      if newval < val
         xold = xcurr; 
         xcurr=newx; 
         val = newval; 
         gradf= newgrad; 
         H = newH;
         Z = [];
         if pnewt == 2, 
            dnewt = newgrad(1:n,2); 
         end
         g = project(A,-gradf(:,1),PT,LZ,UZ,pcolZ,PZ);
         gnrm = norm(g,inf);

         %          Extract the Newton direction?
         if pnewt == 2, 
            dnewt = newgrad(1:n,2); 
         end
      end
      it = it+1; 
      vval(it,1) = val;
   end
   if it > maxcount, 
      ex=4; EXITFLAG = 0;
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
end
if showstat >1,
   display1('final',figtr); 
end
if showstat, 
   xplot(it,vval,vgnrm,vflops,vpos,[],vpcg); 
end

HESSIAN = H;
GRAD = g;
FVAL = val;
if computeLambda
   LAMBDA.eqlin = -A'\gradf;
   LAMBDA.ineqlin = []; LAMBDA.upper = []; LAMBDA.lower = [];
else
   LAMBDA = [];
end
OUTPUT.iterations = it; OUTPUT.funcCount = it;
OUTPUT.cgiterations = npcg;
OUTPUT.firstorderopt = gnrm;
OUTPUT.algorithm = 'large-scale: projected trust-region Newton';
x(:) = xcurr;

