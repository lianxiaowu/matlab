%����
%�ǵ����ڵ��ֵinterp��ʹ��%��7�� ��С������� leastsq(����M����fitfun.m)
%��С������� curvfit(����M����fitfun2.m)
%����ʽ��ֵ����� polyfit 
%����ʽ��ֵ����� lagrange (����M����lagrange.m)


clear;close;
%�ǵ����ڵ��ֵinterp��ʹ��
x=[0.1,0.2,0.15,0,-0.2,0.3];
y=[0.95,0.84,0.86,1.06,1.50,0.72];
xi=-0.2:0.01:0.3;
yi=interp(x,y,xi,'spline');
subplot(2,2,1)
plot(x,y,'o',xi,yi,'k')
title('spline');

%��С������� leastsq(����M����fitfun.m)
c=leastsq('fitfun',[1,1])
f='a*exp(b*t)';
a=c(1);b=c(2);t=xi;
yi=eval(f);
subplot(2,2,2);
plot(x,y,'o',xi,yi,'k')
title('leastsq');    

%��С������� curvfit(����M����fitfun2.m)
x=[0.1,0.2,0.15,0,-0.2,0.3];
y=[0.95,0.84,0.86,1.06,1.50,0.72];
c=curvefit('fitfun2',[1,1],x,y) 

%����ʽ��ֵ����� polyfit 
p=polyfit(x,y,2)
yi=polyval(p,xi);subplot(2,2,3);
plot(x,y,'o',xi,yi,'k');
title('polyfit');

%����ʽ��ֵ����� lagrange (����M����lagrange.m)
yi=lagrange(x,y,xi);
subplot(2,2,4);
plot(x,y,'o',xi,yi,'k');
title('lagrange');

