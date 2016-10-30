function compareStat()
%%
%ԭʼ�������������
load('duck');
oriradio=InOutratio(duck_node,duck_element);
%%
load('..\tetgen result\tetgen_duck');
tetp=node;
tett=elem;
tetn=showComStat(tetp,tett);
tetn=tetn/sum(tetn);
%���������ε�����
tettri=getSurf(tett);
tetradio=InOutratio(tetp,tettri);
%%
load('..\vic result\final_duck');
n=showComStat(p,t);
n=n/sum(n);
%���������ε�����
tri=getSurf(t);
radio=InOutratio(p,tri);

%%
maxy=max([tetn;n])*1.1;


figure(10);
%%
subplot(2,3,4);
oriradio=histc(oriradio,0:0.0125:0.5);
bar(0:0.0125:0.5,oriradio,1,'c');
title('����ģ�ͱ���')

subplot(2,3,2);
bar(0:0.0125:1,tetn,1,'c');
title('iso2mesh��������')
axis([-0.1 1.1 0 maxy]);

subplot(2,3,5);
tetradio=histc(tetradio,0:0.0125:0.5);
bar(0:0.0125:0.5,tetradio,1,'c');
title('iso2mesh���ɱ���')

subplot(2,3,3);
bar(0:0.0125:1,n,1,'c');
title('�����㷨��������')
axis([-0.1 1.1 0 maxy]);

subplot(2,3,6);
radio=histc(radio,0:0.0125:0.5);
bar(0:0.0125:0.5,radio,1,'c');
title('�����㷨������')
end


function n=showComStat(p,t)
%SHOWSTA Summary of this function goes here
%p=3r/R
area=triArea(p,t); %����������ĸ�������
vol=triVolume(p,t);
rin=(3*vol)./sum(area,2);

%�����İ뾶
rout=zeros(size(t,1),1);
for i=1:size(t,1)
    point=[p(t(i,1),:);p(t(i,2),:);p(t(i,3),:);p(t(i,4),:)];
    val=fcircumsphere(point);
    rout(i)=sqrt(val(4));   
end
sta=3*rin./rout;
strcat('�����嵥Ԫ��ƽ������Ϊ��',num2str(mean(sta)))
n=histc(sta,0:0.0125:1);
end
%%
%%�����������
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
%������������
%(1/6)��|det(a ? b, b ? c, c ? d)|, 
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
end

function out=InOutratio(p,tri)
%http://iask.sina.com.cn/b/16818518.html
A=sqrt(sum((p(tri(:,1),:)-p(tri(:,2),:)).^2,2));
B=sqrt(sum((p(tri(:,1),:)-p(tri(:,3),:)).^2,2));
C=sqrt(sum((p(tri(:,2),:)-p(tri(:,3),:)).^2,2));
S=(A+B+C)/2;
temp=(S-A).*(S-B);
temp=temp.*(S-C);
out=4*temp./(A.*B.*C);
strcat('���������ε�ƽ������Ϊ��',num2str(2*mean(out)))
end