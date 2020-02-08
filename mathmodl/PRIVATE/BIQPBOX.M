function[nx,nsig,alpha] = biqpbox(s,c,strg,x,y,sigma,l,u,...
   oval,po,normg,DS,mtxmpy,data,H)
%BIQPBOX Bisection reflective line search for sqpbox
%
%   [nx,nsig,alpha] = BIQPBOX(s,c,strg,x,y,sigma,l,u,...
%                                    oval,po,normg,DS,mtxmpy,data,H)
%   returns the new feasible point nx, the corresponding sign vector nsig,
%   and the step size of the unreflected step, alpha.

%   Copyright (c) 1990-98 by The MathWorks, Inc.
%   $Revision: 1.2 $  $Date: 1998/03/21 16:29:09 $

fac = 10^8; kbnd = 14;
kbnd2 = 4;
lsig = 10^(-8); usig = .9;
strgr = []; strgl = [];

% COMPUTE THE SECOND_ORDER TERM OF THE QUADRATIC APPROXIMATION
k = 0; ss = sigma.*s;
w = DS*ss; ww = feval(mtxmpy,w,H,data); w = DS*ww;
term2 = min(0,.5*ss'*w);

% TRY INITIAL STEP
if po > 0
   alpha = 1; right = 1; sbnd = fac*normg;
else
   alpha = 2; right = 2; sbnd = inf;
end
righty = y+alpha*s;
[rightx,sigr] = reflect(righty,u,l);
w = DS*rightx; ww = feval(mtxmpy,w,H,data); rightw = DS*ww;
rightg = rightw + c; val = rightx'*((.5*rightw +c));
strgr = s'*(sigr.*rightg);
if val < oval + lsig*( alpha*strg + (alpha^2)*term2)
   nx = rightx; nsig = sigr; nstrg = strgr;
else
   %    BISECT UNTIL AN ACCEPTABLE POINT IS FOUND (OR UNTIL k REACHES
   %    UPPER BOUND kbnd). AN INVARIANT IN THE MAIN LOOP IS: AN ACCEPTABLE
   %    RANGE FOR alpha IS CONTAINED IN [left,right].
   left = 0; lefty=y; leftx = x; sigl = sigma; strgl = strg;
   for k=1: kbnd 
      mid = (left + right)/2;
      midy = y + mid*s;
      [midx,sigm] = reflect(midy,u,l);
      w = DS*midx; ww = feval(mtxmpy,w,H,data); midw = DS*ww;
      midg = midw + c;  val = midx'*((.5*midw +c));
      strgmid = s'*(sigm.*midg);
      if val < oval + usig*(mid*strg + (mid^2)*term2) 
         left = mid; strgl = strgmid;
      elseif val > oval + lsig*(mid*strg + (mid^2)*term2)
         right = mid; strgr = strgmid;
      else
         nx = midx; nsig = sigm; alpha = mid; nstrg = strgmid;
         break; 
      end
   end
end
if kbnd2 <= 0, return; end
if k == kbnd
   alpha = 0; nx = x; nsig = sigma;
   return
end
% CONTINUE BISECTION TO FIND A BETTER (AND ACCEPTABLE) POINT
if nstrg < 0
   left = alpha; lefty = y + alpha*s; leftval = val; leftx = nx;
   right = alpha + min(.5,sbnd)*alpha;
   righty = y + right*s; 
   [rightx,sigr] = reflect(righty,u,l);
   w = DS*rightx; ww = feval(mtxmpy,w,H,data); rightw = DS*ww;
   rightg = rightw + c; rightval = rightx'*((.5*rightw +c));
   strgr = s'*(sigr.*rightg);
   % BISECT TO THE RIGHT
   for i=1:kbnd2
      mid = (left + right)/2;
      midy = y + mid*s;
      [midx,sigm] = reflect(midy,u,l);
      w = DS*midx; ww = feval(mtxmpy,w,H,data); midw = DS*ww;
      midg = midw + c;  val = midx'*((.5*midw +c));
      strgmid = s'*(sigm.*midg);
      if (strgmid < 0) & (val < leftval)
         left = mid; leftval = val;
      else
         right = mid; rightval = val;
      end
   end
else  % nstrg >= 0
   right = alpha; righty = y + alpha*s; rightval = val; rightx = nx;
   left  = alpha - min(.5,sbnd)*alpha;
   lefty = y + left*s;
   [leftx,sigl] = reflect(lefty,u,l);
   w = DS*leftx; ww = feval(mtxmpy,w,H,data); leftw = DS*ww;
   leftg = leftw + c; leftval = leftx'*((.5*leftw +c));
   strgl = s'*(sigl.*leftg);
   % BISECT TO THE LEFT
   for i=1:kbnd2
      mid = (left + right)/2;
      midy = y + mid*s;
      [midx,sigm] = reflect(midy,u,l);
      w = DS*midx; ww = feval(mtxmpy,w,H,data); midw = DS*ww;
      midg = midw + c;  val = midx'*((.5*midw +c));
      strgmid = s'*(sigm.*midg);
      if (strgmid > 0) & (val < rightval)
         right = mid; rightval = val;
      else
         left = mid; leftval = val;
      end
   end % for i
end % if nstrg < 0

% CHOOSE THE SMALLER OF THE TWO INTERVAL ENDPOINTS, RETURN.
if leftval <= rightval
   alpha = left; nstrg = strgl;
else
   alpha = right; nstrg = strgr;
end
midy = y + alpha*s;
[nx,nsig] = reflect(midy,u,l);


