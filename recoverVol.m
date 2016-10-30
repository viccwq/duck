function [ output_args ] = recoverVol( node,ele)
%RECOVERVOL Summary of this function goes here
%�ָ�������

%%
Amplify=1;    %����Ŵ�

%%
%�������
node=node*Amplify;      %�Ŵ�����
[ node_new,range,temp ] = sizeFun3dOri(node);
node=node_new;
clear node_new;
%��ȡ�߶εĲ�������
[nodeOri,v123]=barFunc(node,ele);
%����������
out1=int8(zeros(ceil(range(1:2))));
outn=repmat(out1,[1,1,ceil(range(3))]);
figure
for i=1:ceil(range(1))
    YZ=findYZ(i-0.5,nodeOri,v123);
    out1=full(sparse(YZ(:,1),YZ(:,2),ones(size(YZ,1),1),ceil(range(2)),ceil(range(3))));
    outn(i,:,:)=out1;
    h=imshow(out1*252./out1);
%     h=pcolor(double(outn(:,:,i)./outn(:,:,i)));
%     colormap(gray);
%     axis equal;
    pause(0.02);
    delete(h);
    
end
return
end
%%
function [nodeOri,v123] = barFunc(node,ele)
%INSERTPOINT Summary of this function goes here
%   Detailed explanation goes here
%�����������߶εĲ�������

%�ҳ����еı�
pair=zeros(0,2);
localpairs=nchoosek(1:3,2); %�������1-2 1-3 1-4 2-3 ...
for ii=1:size(localpairs,1)     %��ȡ���˽ṹ�ı�
    pair=[pair;ele(:,localpairs(ii,:))];
end
pair=unique(sort(pair,2),'rows');
%V1*t=X-X0     0<=t<=1
%V2*t=Y-Y0     0<=t<=1
%V3*t=Z-Z0     0<=t<=1
nodeOri=node(pair(:,1),:);
v123=node(pair(:,2),:)-node(pair(:,1),:);%xyz�����ϵ���������Inf ���߶δ�ֱ�ڶ�Ӧ��������

return
end
%%
%�ҳ���Ӧ��YZ������
%X0+V1*t=X     0<=t<=1
%Y0+V2*t=Y     0<=t<=1
%Z0+V3*t=Z     0<=t<=1
function YZ=findYZ(X,nodeOri,v123)
%��xEnd=Xʱ
t=(X-nodeOri(:,1))./v123(:,1);
index=((t>=0)&(t<=1));
YZ=repmat(t(index),1,2).*v123(index,2:3)+nodeOri(index,2:3);
YZ=round(YZ);
temp=YZ==0;
YZ(temp)=YZ(temp)+1;
return
end

%%
%�ҳ���Ӧ��XY������
%X0+V1*t=X     0<=t<=1
%Y0+V2*t=Y     0<=t<=1
%Z0+V3*t=Z     0<=t<=1
function XY=findXY(Z,nodeOri,v123)
%��xEnd=Xʱ
t=(Z-nodeOri(:,3))./v123(:,3);
index=((t>=0)&(t<=1));
XY=repmat(t(index),1,2).*v123(index,1:2)+nodeOri(index,1:2);
XY=round(XY);
temp=XY==0;
XY(temp)=XY(temp)+1;
return
end