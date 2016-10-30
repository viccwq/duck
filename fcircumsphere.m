function [ output ] = fcircumsphere( p ) 
%**********************************************
%FCIRCUMSPHERE
%返回四面体的外接球球中心坐标和半径的平方;
%%
%3维离散数据四面体快速生成算法研究
%[ output ] = fcircumsphere( p )
%input:     point(1,1:3)=[x,y,z]is the coordinate of a point in one row,there are four points 
%output:    output=[x,y,z,radius^2];球半径的平方
%%

p=p(1:4,:);
%与p(1,:）相邻棱的向量
a= [p(2,:)-p(1,:);
    p(3,:)-p(1,:);
    p(4,:)-p(1,:)];
%与p(1,:）相邻棱的中点
m= [p(2,:)+p(1,:);
    p(3,:)+p(1,:);
    p(4,:)+p(1,:)];
m=m/2;
%ax=b;
b= [dot(a(1,:),m(1,:));
    dot(a(2,:),m(2,:));
    dot(a(3,:),m(3,:))];
if(det(a)==0)
    error('vic——矩阵的行列是为0');
end
x=(a^-1)*b;
r1=sum((x'-p(1,:)).^2);
% r2=sum((x'-p(2,:)).^2);
% r3=sum((x'-p(3,:)).^2);
% r4=sum((x'-p(4,:)).^2);
output=[x',r1];
end

