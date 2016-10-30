%function MyTest(x,y)
function showTetra(node,element,fix)
%%
% [nodePY,ele]=sortPoin(node,element,fix);
%%
[tri1,tri2]=sotrPoint2(node,element,fix);
% bcol=[.8,.9,1];
bcol=[.001,1,1];
icol=[0.8863 0.2667 0.0863];

%%
close all
h=figure;
set(h,'rend','opengl');
hold on;
% tetramesh(ele,node,'edgecolor','k','FaceAlpha',1,'FaceColor',[0,0.5,0.8]);
% plot3(nodePY(:,1),nodePY(:,2),nodePY(:,3),'.r','MarkerSize',20); 

% h=trimesh(tri1,node(:,1),node(:,2),node(:,3),'facecolor',icol,'edgecolor','k','FaceAlpha',1);
h=trimesh(tri2,node(:,1),node(:,2),node(:,3),'facecolor',bcol,'edgecolor','y','FaceAlpha',0.5);
% plot3(p(:,1),p(:,2),p(:,3),'ok','markerfacecolor','r');

xlabel ( '<--- X --->' );
ylabel ( '<--- Y --->' );
zlabel ( '<--- Z --->' );
title ( 'filename' );

view(3)
axis equal off  vis3d;
rotate3d on
grid on
cameramenu;

return
end
%%
%���ڲ��ĵ�
function [nodePY,ele]=sortPoin(node,element,fix)
%ƽ�ƺ�Ĺ̶���
numFix=size(fix,1);
nodePY=node(1:numFix,:);

%��Բ�����ĵ����
farthest=max(sqrt(nodePY(:,1).^2+nodePY(:,2).^2+nodePY(:,3).^2));

%�����̶����������
temp=element<=numFix;
temp2=temp(:,1)|temp(:,2)|temp(:,3)|temp(:,4);
ele=element(temp2,:);

%ѡȡ����x��0��������
sets=1; %x����
gX=(node(ele(:,1),sets)+node(ele(:,2),sets)+node(ele(:,3),sets)+node(ele(:,4),sets))/4;
sets=2;
gY=(node(ele(:,1),sets)+node(ele(:,2),sets)+node(ele(:,3),sets)+node(ele(:,4),sets))/4;
sets=3;
gZ=(node(ele(:,1),sets)+node(ele(:,2),sets)+node(ele(:,3),sets)+node(ele(:,4),sets))/4;

ele=ele(sqrt(gX.^2+gY.^2+gZ.^2)<farthest,:);


sets=1; %x����
gX=(node(ele(:,1),sets)+node(ele(:,2),sets)+node(ele(:,3),sets)+node(ele(:,4),sets))/4;
ele=ele(gX<0,:);
%

return
end
%%
%����fix������������
function [tri1,tri2]=sotrPoint2(node,element,fix)
%ƽ�ƺ�Ĺ̶���
numFix=size(fix,1);

%�����̶����������
ele=element(any(ismember(element,1:numFix),2),:);

%���������̶����������
% temp=(element<=numFix);
% temp=sum(temp,2);
% temp=(temp==3);
% ele=element(temp,:);

%����ɸѡ
% ele=ele(node(ele(:,1),1)>0,:);

%���ظ���������
faces=[ele(:,[1,2,3]);
       ele(:,[1,2,4]);
       ele(:,[1,3,4]);
       ele(:,[2,3,4])];
faces=sort(faces,2);
[foo,ix,jx]=unique(faces,'rows');
vec=histc(jx,1:max(jx));
qx=find(vec==1);
tri=faces(ix(qx),:);

tri1=tri(all(ismember(tri,1:numFix),2),:);
tri2=setdiff(tri,tri1,'rows');
return
end