%�������Ż�
clear;
%��Լ���Ż�fmins(����M�ļ�optfun1.m)
x=[-1.2,1];
[x,options]=fmins('optfun1',x);
x,options(8)

%Լ���Ż�constr (����M�ļ�optfun2.m)
x0=[20,10,16];
options(13)=1;
[x,options]=constr('optfun2',x0,options,...  
    [-inf,10,-inf]',[inf,20,inf]);
x,-options(8)
