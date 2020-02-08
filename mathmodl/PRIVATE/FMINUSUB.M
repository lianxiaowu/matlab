function [x,FVAL,GRADIENT,HESSIAN,EXITFLAG,OUTPUT] = fminusub(funfcn,x,verbosity,options,Fval,Gval,Hval,varargin)
%FMINUSUB  Finds the minimum of a function of several variables.

%   Copyright (c) 1990-98 by The MathWorks, Inc.
%   $Revision: 1.9 $  $Date: 1998/09/15 21:13:36 $
%   Andy Grace 7-9-90.


% ------------Initialization----------------
HESSIAN = [];
GRADIENT = [];
EXITFLAG = 1;  % assume convergence
OUTPUT = [];
iter = 0;

XOUT=x(:);
% numberOfVariables must be the name of this variable:
numberOfVariables=length(XOUT);    
GRAD=zeros(numberOfVariables,1);
GRAD=Gval;
HESS=Hval;
f=Fval;

n = length(XOUT); numberOfVariables = n;
OLDX=XOUT;
MATX=zeros(3,1);
MATL=[f;0;0];
OLDF=f;
[OLDX,OLDF,HESS]=optinit(XOUT,f,verbosity);  % subfunction
FIRSTF=f;
FIRSTGRAD=GRAD;
FIRSTHESS=HESS;
CHG = 1e-7*abs(XOUT)+1e-7*ones(numberOfVariables,1);
SD = zeros(numberOfVariables,1);
diff = zeros(numberOfVariables,1);
PCNT = 0;
how = '';

gradflag =  strcmp(optimget(options,'GradObj'),'on');
tolX = optimget(options,'tolx');
% lineSearchType of 0 means quadcubic (the default)
lineSearchType = strcmp(optimget(options,'linesearchtype'),'cubicpoly');
hessUpdate = optimget(options,'HessUpdate');
switch optimget(options,'HessUpdate')
case 'bfgs'
   hessUpdate = 0;
case 'dfp'
   hessUpdate = 1;
case 'gillmurray'
   hessUpdate = 3;
case 'steepdesc'
   hessUpdate = 2;
end

tolFun = optimget(options,'tolfun');
DiffMinChange = optimget(options,'diffminchange');
DiffMaxChange = optimget(options,'diffmaxchange');
DerivativeCheck = strcmp(optimget(options,'DerivativeCheck'),'on');
maxFunEvals = optimget(options,'maxfunevals');
maxIter = optimget(options,'maxIter');
% In case the defaults were gathered from calling: optimset('fminsearch'):
if ischar(maxFunEvals)
   maxFunEvals = eval(maxFunEvals);
end

numFunEvals = 0;
numGradEvals = 0;

status =-1;

while status ~= 1
   iter = iter + 1;
   % Work Out Gradients
   if ~(gradflag) | DerivativeCheck 
      GRADFD = GRAD;  % set to correct size
      OLDF=f;
      % Finite difference perturbation levels
      % First check perturbation level is not less than search direction.
      f = find(10*abs(CHG)>abs(SD)); 
      CHG(f) = -0.1*SD(f);
      % Ensure within user-defined limits
      CHG = sign(CHG+eps).*min(max(abs(CHG),DiffMinChange),DiffMaxChange);
      for gcnt=1:numberOfVariables
         XOUT(gcnt,1)=XOUT(gcnt)+CHG(gcnt);
         x(:) = XOUT; 
         f = feval(funfcn{3},x,varargin{:});
         GRADFD(gcnt)=(f-OLDF)/(CHG(gcnt));
         if f < OLDF
            OLDF=f;
         else
            XOUT(gcnt)=XOUT(gcnt)-CHG(gcnt);
         end
      end
      % Try to set difference to 1e-8 for next iteration
      % Add eps for machines that can't handle divide by zero.
      CHG = 1e-8./(GRADFD + eps); 
      f = OLDF;
      numFunEvals=numFunEvals+numberOfVariables;
      % Gradient check 
      if DerivativeCheck == 1 & gradflag
         if isa(funfcn{4},'inline')
            graderr(GRADFD, GRAD, formula(funfcn{4}));
         else
            graderr(GRADFD, GRAD,  funfcn{4});
         end
         DerivativeCheck = 0; 
      else
         GRAD = GRADFD;
      end
   else
      x(:)=XOUT; 
   end
   %---------------Initialization of Search Direction------------------
   if status == -1
      SD=-GRAD;
      FIRSTF=f;
      FIRSTGRAD=GRAD; FIRSTHESS = HESS;
      OLDG=GRAD;
      GDOLD=GRAD'*SD;
      % For initial step-size guess assume the minimum is at zero. 
      currentstepsize = max(0.001, min([1,2*abs(f/GDOLD)]));
      if verbosity>1
         disp([sprintf(' %5.0f       %5.0f    %13.6g  %13.6g   ',iter,numFunEvals,f,currentstepsize),sprintf('%12.3g  ',GDOLD)]);
      end
      XOUT=XOUT+currentstepsize*SD;
      status=4; 
      if lineSearchType==0 
         PCNT=1; 
      end
      
   else
      %-------------Direction Update------------------
      gdnew=GRAD'*SD;
      if verbosity>1, 
         num=[sprintf(' %5.0f       %5.0f    %13.6g  %13.6g   ',iter,numFunEvals,f,currentstepsize),sprintf('%12.3g  ',gdnew)];
      end
      if (gdnew>0 & f>FIRSTF)|~finite(f)
         % Case 1: New function is bigger than last and gradient w.r.t. SD -ve
         %   ...interpolate.
         how='inter';
         [stepsize]=cubici1(f,FIRSTF,gdnew,GDOLD,currentstepsize);
         if stepsize<0|isnan(stepsize), 
            stepsize=currentstepsize/2; 
            how='C1f '; 
         end
         if currentstepsize<0.1& hessUpdate ==0 
            if stepsize*norm(SD)<eps
               stepsize=exp(rand(1,1)-1)-0.1;
               how='RANDOM STEPLENGTH';
               status=0;
            else        
               stepsize=stepsize/2;
            end   
         end      
         currentstepsize=stepsize;
         XOUT=OLDX;
      elseif f<FIRSTF
         [newstep,fbest] =cubici3(f,FIRSTF,gdnew,GDOLD,currentstepsize);
         sk=(XOUT-OLDX)'*(GRAD-OLDG);
         if sk>1e-20
            % Case 2: New function less than old fun. and OK for updating HESS
            %         .... update and calculate new direction.
            how='';   
            if gdnew<0
               how='incstep';
               if newstep<currentstepsize,  
                  newstep=2*currentstepsize+1e-5; 
                  how=[how,' IF']; 
               end
               currentstepsize=min([max([2,1.5*currentstepsize]),1+sk+abs(gdnew)+max([0,currentstepsize-1]), (1.2+0.3*(~lineSearchType))*abs(newstep)]);
            else % gdnew>0
               if currentstepsize>0.9
                  how='int_st';
                  currentstepsize=min([1,abs(newstep)]);
               end
            end %if gdnew
            [HESS,SD]=updatehess(XOUT,OLDX,GRAD,OLDG,HESS,hessUpdate);  % subfunction
            gdnew=GRAD'*SD;
            OLDX=XOUT;
            status=4;
            % Save Variables for next update
            FIRSTF=f;
            FIRSTGRAD=GRAD;
            FIRSTHESS = HESS;
            OLDG=GRAD;
            GDOLD=gdnew;
            % If mixed interpolation set PCNT
            if lineSearchType==0, 
               PCNT=1; MATX=zeros(3,1);  
               MATL(1)=f; 
            end
         elseif gdnew>0 %sk<=0 
            % Case 3: No good for updating HESSIAN .. interpolate or halve step length.
            how='inter_st'; 
            if currentstepsize>0.01
               currentstepsize=0.9*newstep;
               XOUT=OLDX;
            end
            if currentstepsize>1, 
               currentstepsize=1; 
            end
         else  
            % Increase step, replace starting point
            currentstepsize=max([min([newstep-currentstepsize,3]),0.5*currentstepsize]);
            how='incst2';
            OLDX=XOUT;
            FIRSTF=f;
            FIRSTGRAD = GRAD;
            FIRSTHESS = HESS;
            OLDG=GRAD;
            GDOLD=GRAD'*SD;
            OLDX=XOUT;
         end % if sk>
         % Case 4: New function bigger than old but gradient in on
         %         ...reduce step length.
      else %gdnew<0 & F>FIRSTF
         if gdnew<0&f>FIRSTF
            how='red_step';  
            if norm(GRAD-OLDG)<1e-10 
               HESS=eye(numberOfVariables);
            end
            if abs(currentstepsize)<eps
               SD=norm(numberOfVariables,1)*(rand(numberOfVariables,1)-0.5);
               currentstepsize=abs(rand(1,1)-0.5)*1e-6;
               how='RANDOM SD';
            else
               currentstepsize=-currentstepsize/2;
            end
            XOUT=OLDX;
         end %gdnew>0    
      end % if (gdnew>0 & F>FIRSTF)|~finite(F)
      XOUT=XOUT+currentstepsize*SD;
      if isinf(verbosity) 
         disp([num,how])
      elseif verbosity>1 
         disp(num)
      end
   end %----------End of Direction Update-------------------
   
   % Check Termination 
   if max(abs(SD)) < 2*tolX 
      if verbosity > 0
         disp(' ')
         disp('Optimization terminated successfully:')
         disp(' Search direction less than 2*options.TolX')     
      end
      status=1; 
   elseif (-GRAD'*SD) < 2*tolFun & (GRAD'*SD) < 0  % descent direction
      if verbosity > 0
         disp(' ')
         disp('Optimization terminated successfully:')
         disp(' Current search direction is a descent direction, and magnitude of ')
         disp(' directional derivative in search direction less than 2*options.TolFun')     
      end
      status=1; 
      
   elseif numFunEvals > maxFunEvals 
      if verbosity > 0
         disp(' ')
         disp('Maximum number of function evaluations exceeded;')
         disp('   increase options.MaxFunEvals')
      end              
      status=1;
      EXITFLAG = 0;
   elseif iter > maxIter 
      if verbosity > 0
         disp(' ')
         disp('Maximum number of iterations exceeded;')
         disp('   increase options.MaxIter')
      end              
      status=1;
      EXITFLAG = 0;
   else % continue iterating
      % Line search using mixed polynomial interpolation and extrapolation.
      if PCNT~=0 
         while PCNT > 0 & numFunEvals <= maxFunEvals
            x(:) = XOUT;
            f = feval(funfcn{3},x,varargin{:});
            numFunEvals=numFunEvals+1;
            [PCNT,MATL,MATX,steplen,f, how]=searchq(PCNT,f,OLDX,MATL,MATX,SD,GDOLD,currentstepsize, how);
            currentstepsize=steplen;
            XOUT=OLDX+steplen*SD;
         end % end while
         if numFunEvals>maxFunEvals 
            if verbosity > 0
               disp(' ')
               disp('Maximum number of function evaluations exceeded;')
               disp('   increase options.MaxFunEvals')
            end              
            status=1; 
            EXITFLAG = 0;
         end
      end % if PCNT~=0
   end % if max...
   
   x(:)=XOUT;
   switch funfcn{1}
   case 'fun'
      f = feval(funfcn{3},x,varargin{:});
   case 'fungrad'
      [f,GRAD(:)] = feval(funfcn{3},x,varargin{:});
      numGradEvals=numGradEvals+1;
   case 'fun_then_grad'
      f = feval(funfcn{3},x,varargin{:}); 
      GRAD(:) = feval(funfcn{4},x,varargin{:});
      numGradEvals=numGradEvals+1;
   otherwise
      error('Undefined calltype in FMINUNC');
   end
   numFunEvals=numFunEvals+1;
end  % while

if f > FIRSTF
     FVAL = FIRSTF;
   GRADIENT = FIRSTGRAD;
   HESSIAN = FIRSTHESS;
   x(:)=OLDX;
else
   FVAL = f;
   GRADIENT = GRAD;
   HESSIAN = HESS;
end

OUTPUT.iterations = iter;
OUTPUT.funcCount = numFunEvals;
OUTPUT.stepsize=currentstepsize;
OUTPUT.firstorderopt = norm(GRAD,inf);
OUTPUT.algorithm='medium-scale: Quasi-Newton line search';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [xold,fold,invhess]=optinit(xnew,fnew,verbosity)
%OPTINT Function to initialize FMINU routine.

%   Copyright (c) 1990-98 by The MathWorks, Inc.
%   $Revision: 1.9 $  $Date: 1998/09/15 21:13:36 $
%   Andy Grace 7-9-90.
lenx=length(xnew);
invhess=eye(lenx);  
xold=xnew;
fold=fnew;
if verbosity > 1
   if isinf(verbosity)
      header = sprintf(['\n                                                      Directional \n',...
            ' Iteration  Func-count     f(x)         Step-size      derivative  Procedure ']);
      
   else
      header = sprintf(['\n                                                      Directional \n',...
            ' Iteration  Func-count     f(x)         Step-size      derivative ']);
   end
   disp(header)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [invhess,directn]=updatehess(xnew,xold,gradxnew,gradxold,invhess,hessUpdate)
%UPDHESS Performs the Inverse Hessian Update.
%   Returns direction of search for use with 
%   unconstrained optimization problems. 

%   Copyright (c) 1990-98 by The MathWorks, Inc.
%   $Revision: 1.9 $  $Date: 1998/09/15 21:13:36 $
%   Andy Grace 7-9-90.

u=xnew-xold;
v=gradxnew-gradxold;
if hessUpdate==0
% The BFGS Hessian Update formula:
    invhess=invhess + v*v'/(v'*u)  -invhess*u*u'*invhess'/(u'*invhess*u);
    directn=-invhess\gradxnew;

elseif hessUpdate==1

% The DFP formula
    a=u*u'/(u'*v);
    b=-invhess*v*v'*invhess'/(v'*invhess*v);
    invhess=invhess + a + b;
    directn=-invhess*gradxnew;

elseif hessUpdate==3

% A formula given by Gill and Murray
    a = 1/(v'*u);
    invhess=invhess - a*(invhess*v*u'+u*v'*invhess)+a*(1+v'*invhess*v*a)*u*u' ;
    directn=-invhess*gradxnew;
elseif hessUpdate==2
% Steepest Descent
    directn=-gradxnew;
end
       
