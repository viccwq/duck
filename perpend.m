function [ dist,vec ] = perpend( p,bar1,bar2 )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%%
%�󹫴���
%input p    ����Ľڵ� m*3
%input p0   �߶�0�����
%input p0e  �߶�0���յ�
%input p1   �߶�0�����
%input p1e  �߶�0���յ�

%%
p0=p(bar1(:,1),:);
p0e=p(bar1(:,2),:);
p1=p(bar2(:,1),:);
p1e=p(bar2(:,2),:);
%����ֱ�ߵ�����
d0=p0e-p0;
d1=p1e-p1;
%
a=dot(d0,d0,2);
b=-1*dot(d0,d1,2);
c=dot(d1,d1,2);
d=dot(d0,p0-p1,2);
e=dot(d0,p0-p1,2);
%����
sc=(b.*e-c.*d)./(a.*c-b.*b);
tc=(b.*d-a.*e)./(a.*c-b.*b);
Q0=p0+repmat(sc,1,3).*d0;
Q1=p1+repmat(tc,1,3).*d1;
%���س���
dist=sqrt(sum((Q1-Q0).^2,2));
%���ص�λ����
vec=Q0-Q1;
vec=vec./repmat(dist,1,3);



end

