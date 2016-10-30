function kill_surf_kite(stdn)
load pANDt3d
StaQShort(p,t,stdn);
end

function StaQShort(p,t,stdn)
%SHOWSTA Summary of this function goes here
%p=3r/R
area=triArea(p,t); %求的四面体四个面的面积
vol=triVolume(p,t);
rin=(3*vol)./sum(area,2);

%外接球的半径
rout=zeros(size(t,1),1);
for i=1:size(t,1)
    point=[p(t(i,1),:);p(t(i,2),:);p(t(i,3),:);p(t(i,4),:)];
    val=fcircumsphere(point);
    rout(i)=sqrt(val(4));   
end
sta=3*rin./rout;
figure(4);
subplot(2,2,1); 
n=histc(sta,0:0.0125:1);
bar(0:0.0125:1,n,1,'c');
title('未删除表面kite后的结果');

% n=histc(sta,min(sta):0.01:max(sta));
% bar(min(sta):0.01:max(sta),n,1,'c');


%% 显示面
subplot(2,2,2); 
tri=getSurf(t);
tri1=tri(p(tri(:,1),3)>40,:);
tri1=sort(tri1,2);
trimesh(tri1,p(:,1),p(:,2),p(:,3),'facecolor',[0,0.5,0.8],'edgecolor','y','FaceAlpha',0.0);
title('not kill deteriorated kite');
%% 显示四边形
% stdn=0.2;
t_kite=t(sta<stdn,:);
tri2=sort(t_kite(:,[1,2,3]),2);
index=ismember(tri2, tri,'rows');
tri2=sort(t_kite(:,[1,2,4]),2);
index=ismember(tri2, tri,'rows')|index;
tri2=sort(t_kite(:,[1,3,4]),2);
index=ismember(tri2, tri,'rows')|index;
tri2=sort(t_kite(:,[2,3,4]),2);
index=ismember(tri2, tri,'rows')|index;
hold on
tetramesh(t_kite(index,:),p,'edgecolor','k','FaceAlpha',1,'FaceColor','r');
tetramesh(t_kite(~index,:),p,'edgecolor','k','FaceAlpha',1,'FaceColor','c');
inner=sum(~index);
string=['比值小于',num2str(stdn),'的单元有',num2str(inner),'/',num2str(sum(sta<stdn)),'个'];
% title(strcat('内、外接球的半径比，',string));
disp(string);
view(3)
axis equal off  vis3d;
rotate3d on
cameramenu
% disp('press any key to continue');

%% 删除表面的kite
t1=t(sta>=stdn,:);%sta>=stdn的四面体
t2=t_kite(~index,:);%取出表面的kite后所剩的四面体

sta1=sta(sta>=stdn);
sta2=sta(sta<stdn);
sta2=sta2(~index,:);

subplot(2,2,3);  
n=histc([sta1;sta2],0:0.0125:1);
bar(0:0.0125:1,n,1,'c');
title('删除表面kite后的结果');

subplot(2,2,4); 
tri=getSurf([t1;t2]);
tri1=tri(p(tri(:,1),3)>40,:);
tri1=sort(tri1,2);
trimesh(tri1,p(:,1),p(:,2),p(:,3),'facecolor','c','edgecolor','k','FaceAlpha',1);
hold on
tetramesh(t_kite(~index,:),p,'edgecolor','k','FaceAlpha',1,'FaceColor','c');
view(3)
axis equal off  vis3d;
rotate3d on
cameramenu
end
%%
%%求四面体面积
function [Tarea]=triArea(p,t)
faces=[t(:,[1,2,3]);
       t(:,[1,2,4]);
       t(:,[1,3,4]);
       t(:,[2,3,4])];
faces=sort(faces,2);
[facesu,ix,jx]=unique(faces,'rows');
vector1=p(facesu(:,2),:)-p(facesu(:,1),:);
vector2=p(facesu(:,3),:)-p(facesu(:,1),:);
area=cross(vector1,vector2,2);
area=sqrt(sum(area.^2,2))/2;
area=area(jx);
Tarea=reshape(area,size(t,1),[]);
end
%%
%求四面体的体积
%(1/6)·|det(a ? b, b ? c, c ? d)|, 
function vol=triVolume(p,t)
v1=p(t(:,2),:)-p(t(:,1),:);
v2=p(t(:,3),:)-p(t(:,1),:);
v3=p(t(:,4),:)-p(t(:,1),:);
tSize=size(t,1);
vol=zeros(tSize,1);
for i=1:tSize
    vol(i)=abs(det([v1(i,:);v2(i,:);v3(i,:)]))/6;
end
end


function tri=getSurf(ele)
%求不重复的三角形
faces=[ele(:,[1,2,3]);
       ele(:,[1,2,4]);
       ele(:,[1,3,4]);
       ele(:,[2,3,4])];
faces=sort(faces,2);
[foo,ix,jx]=unique(faces,'rows');
vec=histc(jx,1:max(jx));
qx=find(vec==1);
tri=faces(ix(qx),:);
end
