function msg = diagnose(caller,OUTPUT,gradflag,hessflag,constflag,gradconstflag,line_search,OPTIONS,XOUT,non_eq,...
   non_ineq,lin_eq,lin_ineq,LB,UB,funfcn,confcn,f,GRAD,HESS,c,ceq,cGRAD,ceqGRAD);
% DIAGNOSE prints diagnostic information about the function to be minimized
%    or solved.

%   Copyright (c) 1990-98 by The MathWorks, Inc.
%   $Revision: 1.8 $ $Date: 1998/09/30 18:17:18 $

msg = [];

pstr = sprintf('\n%s\n%s\n',...
   '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%',...
   '   Diagnostic Information ');
disp(pstr)

if ~isempty(funfcn{1})
   if isa(funfcn{3},'inline')
      funformula = formula(funfcn{3});
   else 
      funformula = funfcn{3};
   end
   if isa(funfcn{4},'inline')
      gradformula = formula(funfcn{4});
   else 
      gradformula = funfcn{4};
   end
   
   if isa(funfcn{5},'inline')
      hessformula = formula(funfcn{5});
   else 
      hessformula = funfcn{5};
   end
end

if ~isempty(confcn{1})
   if isa(confcn{3},'inline')
      conformula = formula(confcn{3});
   else 
      conformula = confcn{3};
   end
   if isa(confcn{4},'inline')
      gradcformula = formula(confcn{4});
   else 
      gradcformula = confcn{4};
   end
else
   conformula = '';
   gradcformula = '';
end

disp(['Number of variables: ', int2str(length(XOUT)),sprintf('\n')])
if ~isempty(funfcn{1})
   disp('Functions ')
   switch funfcn{1}
   case 'fun'
      % display 
      disp([' Objective:                            ',funformula]);
      
   case 'fungrad'
      if gradflag
         disp([' Objective and gradient:               ',funformula]);
      else
         disp([' Objective:                            ',funformula]);
         disp( '   (set OPTIONS.GradObj=''on'' to use user provided gradient function)') 
      end
      
   case 'fungradhess'
      if gradflag & hessflag
         disp([' Objective, gradient and Hessian:      ',funformula]);
      elseif gradflag
         disp([' Objective and gradient:               ',funformula]);
         disp( '   (set OPTIONS.Hessian to ''on'' to use user provided Hessian function)') 
      else
         disp([' Objective:                            ',funformula]);
         disp( '   (set OPTIONS.GradObj=''on'' to use user provided gradient function)')
         disp( '   (set OPTIONS.Hessian to ''on'' to use user provided Hessian function)') 
      end
      
      
   case 'fun_then_grad'
      disp([' Objective:                            ',funformula]);
      if gradflag
         disp([' Gradient:                             ',gradformula]);
      end
      if hessflag
         disp('-->Ignoring OPTIONS.Hessian --no user Hessian function provided')
      end
      
   case 'fun_then_grad_then_hess'
      disp([' Objective:                            ',funformula]);
      if gradflag & hessflag
         disp([' Gradient:                             ',gradformula]);
         disp([' Hessian:                              ',hessformula]);
      elseif gradflag
         disp([' Gradient:                             ',gradformula]);
      end   
   otherwise
      
   end
   
   if ~gradflag
      disp(' Gradient:                             finite-differencing')
   end
   % shape of grad
   
   if ~hessflag & (isequal('fmincon',caller) | isequal('fminunc',caller))
      disp(' Hessian:                              finite-differencing (or Quasi-Newton)')
   end
   % shape of hess
   
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~isempty(confcn{1})
   switch confcn{1}
      
   case 'fun'
      disp([' Nonlinear constraints:                ',conformula]);
   case 'fungrad'
      if gradconstflag
         disp([' Nonlinear constraints and gradient:   ',conformula]);
      else
         disp([' Nonlinear constraints:                ',conformula]);
         disp( '   (set OPTIONS.GradConstr to ''on'' to use user provided gradient of constraints function)') 
      end
      
   case 'fun_then_grad'
      disp([' Nonlinear constraints:                ',conformula]);
      if gradconstflag
         disp([' Nonlinear constraints gradient:       ',gradcformula]);
      end
      
   otherwise
      
   end
   
   if ~constflag
      disp(' Nonlinear constraints:                finite-differencing')
   end
   if ~gradconstflag
      
      disp(' Gradient of nonlinear constraints:    finite-differencing')
   end
   disp([sprintf('\n'),'Constraints'])  
   disp([' Number of nonlinear inequality constraints: ',int2str(non_ineq)])
   disp([' Number of nonlinear equality constraints:   ',int2str(non_eq)])
   
elseif isequal(caller,'fmincon') | isequal(caller,'fminimax') | ...
      isequal(caller,'fgoalattain') | isequal(caller,'fseminf')
   disp([sprintf('\n'),'Constraints'])
   disp(' Nonlinear constraints:             do not exist')
   
end

disp(' ')


switch caller
   
case {'fmincon','linprog','quadprog','lsqlin','fminimax','fseminf','fgoalattain'}
   disp([' Number of linear inequality constraints:    ',int2str(lin_ineq)])
   disp([' Number of linear equality constraints:      ',int2str(lin_eq)])
   disp([' Number of lower bound constraints:          ',int2str(nnz(~isinf(LB)))])
   disp([' Number of upper bound constraints:          ',int2str(nnz(~isinf(UB)))])
case {'lsqcurvefit','lsqnonlin'}
   disp([' Number of lower bound constraints:          ',int2str(nnz(~isinf(LB)))])
   disp([' Number of upper bound constraints:          ',int2str(nnz(~isinf(UB)))])
case {'fsolve','fminunc','fsolves'}
otherwise
end

if ~isempty(OUTPUT)
temp = sprintf('\n%s\n   %s\n','Algorithm selected',OUTPUT.algorithm);
disp(temp)
end

pstr = sprintf('\n%s\n%s\n',...
   '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%',...
   ' End diagnostic information ');
disp(pstr)
