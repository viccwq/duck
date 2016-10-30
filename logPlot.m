
x=1:0.2:40;
a=1.2;
hold on
plot(x,log(x)/log(a),'r');
plot(x,1./(x*log(a)),'b');
plot(x,1,'k');

X=1/log(a);
plot(X,-5:0.2:15,'k');
plot(x,log(X)/log(a),'k');

x=-1:0.2:20;
plot(x,a.^x,'r');
plot(x,(a-0.05).^x,'c');
plot(x,(a-0.02).^x,'y');
axis equal
