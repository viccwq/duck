function [ dist,vec ] = perpend( p,bar1,bar2 )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%%
%求公垂线
%input p    输入的节点 m*3
%input p0   线段0的起点
%input p0e  线段0的终点
%input p1   线段0的起点
%input p1e  线段0的终点

%%
p0=p(bar1(:,1),:);
p0e=p(bar1(:,2),:);
p1=p(bar2(:,1),:);
p1e=p(bar2(:,2),:);
%两条直线的向量
d0=p0e-p0;
d1=p1e-p1;
%
a=dot(d0,d0,2);
b=-1*dot(d0,d1,2);
c=dot(d1,d1,2);
d=dot(d0,p0-p1,2);
e=dot(d0,p0-p1,2);
%参数
sc=(b.*e-c.*d)./(a.*c-b.*b);
tc=(b.*d-a.*e)./(a.*c-b.*b);
Q0=p0+repmat(sc,1,3).*d0;
Q1=p1+repmat(tc,1,3).*d1;
%返回长度
dist=sqrt(sum((Q1-Q0).^2,2));
%返回单位方向
vec=Q0-Q1;
vec=vec./repmat(dist,1,3);



end

