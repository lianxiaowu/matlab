function prim(pp)
%PRIM (����С������)
%����С�������㷨,ͨ��prim�㷨��������,��������Ӧͼ��. 
%�÷�:
%	�����������: 
%		map=[���1 �յ�1 �߳�1;���2 �յ�2  �߳�2;............;���n �յ�n �߳�n]
%	����[out,len]=kruskal(map)��������
%����˵��
%	map----3���ڽӾ���,ÿ�б�ʾһ����.��һ�б�ʾ���,�ڶ��б�ʾ�յ�,�����б�ʾ�߳�
%	out---�������:[��� �յ�]
%	len---������������ܳ���
%
%����  
%		clear;map=[1 2 30;2 4 5;3 2 6;4 3 1;1 4 20;1 5 30]
%		[out,len]=prim(map)
%
%���㷨������VC++6.0����prim1.c���ɵ�MEX�ļ�prim1.dll�����С������,
%	��ʾ��������ֵ�Ͻ�(Ĭ��10000)
%
%See also DIJKSTRA,LPINT,DP,BNBGUI,BNB18,

%By S. J. Lu 2000

mx=0;
[m n]=size(pp);
for i=1:m
   if (pp(i,1)>mx),
      mx=pp(i,1);
   end
   if (pp(i,2)>mx),
      mx=pp(i,2);
   end
end   
for i=1:mx
   for j=1:mx
      op(i,j)=10000;
   end   
end   
for i=1:m
   op(pp(i,1),pp(i,2))=pp(i,3);
   op(pp(i,2),pp(i,1))=pp(i,3);
end
[out,len]=prim1(op);
out=out';

figure;
 set(gcf,'numbertitle','off');
 set(gcf,'name','Prim');
 set(gca,'visible','off');
 axis square;
 hold on;
 
 b=linspace(0,2*pi,mx+1);
 b1=10*sin(b);
 b2=10*cos(b);
 plot(b1,b2,'ko');
 
 hh=char(49:48+mx);
 for i=1:mx
    text(b1(i)+1,b2(i),hh(i));
 end;
 for j=1:m
   for i=1:2
    c1(i)=b1(pp(j,i));
    c2(i)=b2(pp(j,i));
    line(c1,c2);   
    end;
 end;

 d1=0;d2=0;
 k1=length(out);
 for i=1:k1
  for j=1:2
      d1(j)=b1(out(i,j));
      d2(j)=b2(out(i,j));
  end
    h=plot(d1,d2,'r');
    set(h,'linewidth',2);
 end 
 legend(h,'����Ϊ������');