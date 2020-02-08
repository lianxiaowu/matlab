function [x,OPTIONS,CostFunction,JACOB] = nlsqold(FUN,x,OPTIONS,GRADFUN,varargin)
%NLSQ Solves non-linear least squares problems.
%   NLSQ is the core code for solving problems of the form:
%   min  sum {FUN(X).^2}    where FUN and X may be vectors or matrices.   
%             x
%
%   [X,OPTIONS,F,J]=NLSQ('FUN',X0,OPTIONS,'GRADFUN',P1,P2,..) 
%   starts at the matrix X0 and finds a minimum to the
%   sum of squares of the functions described in FUN. FUN is usually
%   an M-file which returns a vector of objective functions: F=FUN(X).
%   NOTE: FUN should return FUN(X) and not the sum-of-squares 
%   sum(FUN(X).^2)) (FUN(X) is summed and squared implicitly in
%   the algorithm).
%
%   OPTIONS is a vector of optional parameters to be defined. 
%   OPTIONS(2) is a measure of the precision required for the 
%   values of X at the solution. OPTIONS(3) is a measure of the precision
%   required of the objective function at the solution. See HELP FOPTIONS. 
%
%   GRADFUN is a function to be entered which returns the partial 
%   derivatives of the functions, dF/dX, (stored in columns) at 
%   the point X: gf = GRADFUN(X).
%
%   P1, P2, ..., are problem-dependent parameters passed directly to 
%   the functions FUN and GRADFUN: FUN(X,P1,P2,...) and 
%   GRADFUN(X,P1,P2,...).  Pass empty matrices for OPTIONS and 'GRADFUN' 
%   to use the default values.
%
%   F is the value of FUN(X) at the solution X, and J the Jacobian of 
%   the function FUN at the solution. 

%   Copyright (c) 1990-98 by The MathWorks, Inc.
%   $Revision: 1.1 $  $Date: 1998/03/03 00:24:22 $
%   Andy Grace 7-9-90.

%   The default algorithm is the Levenberg-Marquardt method with a 
%   mixed quadratic and cubic line search procedure.  A Gauss-Newton
%   method is selected by setting  OPTIONS(5)=1. 
%

% ------------Initialization----------------


XOUT = x(:);
[nvars] = length(XOUT);
how = [];

% Global parameters for outside control of leastsq
% OPT_STOP is used for prematurely stopping the optimization
% OPT_STEP is set to 1 during major (non-gradient finding) iterations
%          set to 0 during gradient finding and 2 during line search
%          this can be useful for plotting etc.
global OPT_STOP OPT_STEP;
OPT_STEP = 1;
OPT_STOP = 0;


CostFunction = feval(FUN,x,varargin{:});
CostFunction = CostFunction(:);
OPT_STEP = 0;  % No longer a major step

nfun=length(CostFunction);
GRAD=zeros(length(XOUT),nfun);
MATX=zeros(3,1);
MATL=[CostFunction'*CostFunction;0;0];
FIRSTF=CostFunction'*CostFunction;
[OLDX,OLDF,OPTIONS]=lsint(XOUT,CostFunction,OPTIONS);
PCNT = 0;
EstSum=0.5;
GradFactor=0; 
CHG = 1e-7*abs(XOUT)+1e-7*ones(nvars,1);

OPTIONS(10)=1;
status=-1;


while status~=1  & OPT_STOP == 0
   
   % Work Out Gradients
   if isempty(GRADFUN) | OPTIONS(9)
      OLDF=CostFunction;
      CHG = sign(CHG+eps).*min(max(abs(CHG),OPTIONS(16)),OPTIONS(17));
      for gcnt=1:nvars
         temp = XOUT(gcnt);
         XOUT(gcnt) = temp +CHG(gcnt);
         x(:) = XOUT;
         CostFunction(:) = feval(FUN,x,varargin{:});
         OPT_STEP = 0; % We're in gradient finding mode
         GRAD(gcnt,:)=(CostFunction-OLDF)'/(CHG(gcnt));
         XOUT(gcnt) = temp;
      end
      CostFunction = OLDF;
      OPTIONS(10)=OPTIONS(10)+nvars;
      % Gradient check
      if OPTIONS(9) == 1 & ~isempty(GRADFUN)
         GRADFD = GRAD;
         x(:)=XOUT; 
         GRAD = feval(GRADFUN,x,varargin{:});
         if isa(GRADFUN,'inline')
            graderr(GRADFD, GRAD, formula(GRADFUN));
         else
            graderr(GRADFD, GRAD,  GRADFUN);
         end
         
         OPTIONS(9) = 0;
      end
   else
      x(:) = XOUT;
      OPTIONS(11)=OPTIONS(11)+1;
      GRAD = feval(GRADFUN,x,varargin{:});
   end
   % Try to set difference to 1e-8 for next iteration
   if nfun==1
      if isequal(GRAD,0)
         CHG = Inf;
      else
         CHG = nfun*1e-8./GRAD;
      end
   else
      sumabsGRAD = sum(abs(GRAD)')';
      ii = (sumabsGRAD == 0);
      CHG(ii) = Inf;
      CHG(~ii) = nfun*1e-8./sumabsGRAD(~ii);
   end
   
   OPT_STEP = 2; % Entering line search    
   
   GradF = 2*GRAD*CostFunction;
   NewF = CostFunction'*CostFunction;
   %---------------Initialization of Search Direction------------------
   if status==-1
      if cond(GRAD)>1e8
         SD=-(GRAD*GRAD'+(norm(GRAD)+1)*(eye(nvars,nvars)))\(GRAD*CostFunction);
         if OPTIONS(5)==0, 
            GradFactor=norm(GRAD)+1; 
         end
         how='COND';
      else
         % SD=GRAD'\(GRAD'*X-F)-X;
         SD=-(GRAD*GRAD'+GradFactor*(eye(nvars,nvars)))\(GRAD*CostFunction);
      end
      FIRSTF=NewF;
      OLDG=GRAD;
      GDOLD=GradF'*SD;
      % OPTIONS(18) controls the initial starting step-size.
      % If OPTIONS(18) has been set externally then it will
      % be non-zero, otherwise set to 1.
      if OPTIONS(18) == 0, 
         OPTIONS(18)=1; 
      end
      if OPTIONS(1)>0
         disp([sprintf('%5.0f %12.6g %12.3g ',OPTIONS(10),NewF,OPTIONS(18)),sprintf('%12.3g  ',GDOLD)]);
      end
      XOUT=XOUT+OPTIONS(18)*SD;
      if OPTIONS(5)==0
         newf=GRAD'*SD+CostFunction;
         GradFactor=newf'*newf;
         SD=-(GRAD*GRAD'+GradFactor*(eye(nvars,nvars)))\(GRAD*CostFunction); 
      end
      newf=GRAD'*SD+CostFunction;
      XOUT=XOUT+OPTIONS(18)*SD;
      EstSum=newf'*newf;
      status=0;
      if OPTIONS(7)==0; 
         PCNT=1; 
      end
      
   else
      %-------------Direction Update------------------
      gdnew=GradF'*SD;
      if OPTIONS(1)>0, 
         num=[sprintf('%5.0f %12.6g %12.3g ',OPTIONS(10),NewF,OPTIONS(18)),sprintf('%12.3g  ',gdnew)];
      end
      if gdnew>0 & NewF>FIRSTF
         % Case 1: New function is bigger than last and gradient w.r.t. SD -ve
         % ... interpolate. 
         how='inter';
         [stepsize]=cubici1(NewF,FIRSTF,gdnew,GDOLD,OPTIONS(18));
         OPTIONS(18)=0.9*stepsize;
      elseif NewF<FIRSTF
         %  New function less than old fun. and OK for updating 
         %         .... update and calculate new direction. 
         [newstep,fbest] =cubici3(NewF,FIRSTF,gdnew,GDOLD,OPTIONS(18));
         if fbest>NewF,
            fbest=0.9*NewF; 
         end 
         if gdnew<0
            how='incstep';
            if newstep<OPTIONS(18),  
               newstep=(2*OPTIONS(18)+1e-4); how=[how,'IF']; 
            end
            OPTIONS(18)=abs(newstep);
         else 
            if OPTIONS(18)>0.9
               how='int_step';
               OPTIONS(18)=min([1,abs(newstep)]);
            end
         end
         % SET DIRECTION.      
         % Gauss-Newton Method    
         temp=1;
         if OPTIONS(5)==1 
            if OPTIONS(18)>1e-8 & cond(GRAD)<1e8
               SD=GRAD'\(GRAD'*XOUT-CostFunction)-XOUT;
               if SD'*GradF>eps,
                  how='ERROR- GN not descent direction',  
               end
               temp=0;
            else
               if OPTIONS(1) > 0
                  disp('Conditioning of Gradient Poor - Switching To LM method')
               end
               how='CHG2LM';
               OPTIONS(5)=0;
               OPTIONS(18)=abs(OPTIONS(18));               
            end
         end
         
         if (temp)      
            % Levenberg_marquardt Method N.B. EstSum is the estimated sum of squares.
            %                                 GradFactor is the value of lambda.
            % Estimated Residual:
            if EstSum>fbest
               GradFactor=GradFactor/(1+OPTIONS(18));
            else
               GradFactor=GradFactor+(fbest-EstSum)/(OPTIONS(18)+eps);
            end
            SD=-(GRAD*GRAD'+GradFactor*(eye(nvars,nvars)))\(GRAD*CostFunction); 
            OPTIONS(18)=1; 
            estf=GRAD'*SD+CostFunction;
            EstSum=estf'*estf;
            if OPTIONS(1)>0, 
               num=[num,sprintf('%12.6g ',GradFactor)]; 
            end
         end
         gdnew=GradF'*SD;
         
         OLDX=XOUT;
         % Save Variables
         FIRSTF=NewF;
         OLDG=GradF;
         GDOLD=gdnew;    
         
         % If quadratic interpolation set PCNT
         if OPTIONS(7)==0, 
            PCNT=1; MATX=zeros(3,1); MATL(1)=NewF; 
         end
      else 
         % Halve Step-length
         how='Red_Step';
         if NewF==FIRSTF,
            if OPTIONS(1)>0,
               disp('No improvement in search direction: Terminating')
            end
            status=1;
         else
            OPTIONS(18)=OPTIONS(18)/8;
            if OPTIONS(18)<1e-8
               OPTIONS(18)=-OPTIONS(18);
            end
         end
      end
      XOUT=OLDX+OPTIONS(18)*SD;
      if isinf(OPTIONS(1))
         disp([num,how])
      elseif OPTIONS(1)>0
         disp(num)
      end
      
   end %----------End of Direction Update-------------------
   if OPTIONS(7)==0, 
      PCNT=1; MATX=zeros(3,1);  MATL(1)=NewF; 
   end
   % Check Termination 
   if max(abs(SD))< OPTIONS(2) & (GradF'*SD) < OPTIONS(3) & max(abs(GradF)) < 10*(OPTIONS(3)+OPTIONS(2))
      if OPTIONS(1) > 0
         disp('Optimization Terminated Successfully')  
      end
      status=1; 
   elseif OPTIONS(10)>OPTIONS(14)
      disp('maximum number of iterations has been exceeded');
      if OPTIONS(1)>0
         disp('Increase OPTIONS(14)')
      end
      status=1;
   else
      
      % Line search using mixed polynomial interpolation and extrapolation.
      if PCNT~=0
         % initialize OX and OLDF 
         OX = XOUT; OLDF = CostFunction;
         while PCNT > 0 & OPT_STOP == 0
            x(:) = XOUT; 
            CostFunction(:) = feval(FUN,x,varargin{:});
            OPTIONS(10)=OPTIONS(10)+1;
            NewF = CostFunction'*CostFunction;
            % <= used in case when no improvement found.
            if NewF <= OLDF'*OLDF, 
               OX = XOUT; 
               OLDF=CostFunction; 
            end
            [PCNT,MATL,MATX,steplen,NewF,how]=searchq(PCNT,NewF,OLDX,MATL,MATX,SD,GDOLD,OPTIONS(18),how);
            OPTIONS(18)=steplen;
            XOUT=OLDX+steplen*SD;
            if NewF==FIRSTF,  
               PCNT=0; 
            end
         end
         XOUT = OX;
         CostFunction=OLDF;
      else
         x(:)=XOUT; 
         CostFunction(:) = feval(FUN,x,varargin{:});
         OPTIONS(10)=OPTIONS(10)+1;
      end
   end
   OPT_STEP = 1; % Call the next iteration a major step for plotting etc.
end
OPTIONS(8) = NewF;
XOUT=OLDX;
x(:)=XOUT;
JACOB = GRAD.';
%--end of leastsq--



