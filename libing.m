clear all;clc;
n=7;
hold on;
t=linspace(0,2*pi,7);      
X=cos(t);
Y=sin(t);
Y(abs(Y)<1E-10)=0;
plot(X,Y,'r');axis equal;
axis([-1.2 1.2 -1.2 1.2]); 
k1=(Y(1)-Y(2))/(X(1)-X(2));
k3=-k1;
k4=k1;
k6=-k1;
%��������
point=rand(100,2)*2-1;
x=point(:,1);y=point(:,2);
y1=(x-1)*k1-y;
y2=Y(2)-y;
y3=(x+1)*k3-y;
y4=(x+1)*k4-y;
y5=Y(5)-y;
y6=(x-1)*k6-y;
y1=y1>0;
y2=y2>0;
y3=y3>0;
y4=y4<0;
y5=y5<0;
y6=y6<0;
mark=y1&y2&y3&y4&y5&y6;
%Բ��
r=0.2; %Բ�İ뾶0~1
t=linspace(0,2*pi,100);      
X=cos(t)*r;
Y=sin(t)*r;
Y(abs(Y)<1E-10)=0;
plot(X,Y,'b');
mark=mark&(sum(point.^2,2)>r^2);
x=x(mark);
y=y(mark);
%ȡʮ����
x(11:end)=[];
y(11:end)=[];
plot(x,y,'.k','MarkerSize',5);
