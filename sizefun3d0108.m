function [ x,y,vic_h2 ] = sizefun3d( )
%SIZEFUN Summary of this function goes here
%   Detailed explanation goes here
%%
close all
load cow_node;
load cow_element;
node=cow_node;
element=cow_element;
clear cow_node cow_element;

Desc=1.2;     %笛卡尔网格大小 Descartes
DDensity=2;%笛卡尔网格密度，为图形最小尺寸的百分比
Alpha=20;      %影响系数矩阵 底数
Beta=12;        %计算尺寸函数 底数
Amplify=10;    %坐标放大

node=node*Amplify;   %扩大坐标
eledist=elementdis(node,element);
% eledist=eledist/max(eledist);
%%
%获取点的范围,将图形的重心平移到原点
range=(minmax(node'))';
node=node-repmat((range(1,:)+range(2,:))/2,size(node,1),1);
range=range-repmat((range(1,:)+range(2,:))/2,2,1);

%笛卡尔网格
[x,y]=meshgrid( (-range(2,1)*Desc):(min(eledist)*DDensity):(range(2,1)*Desc),...
                (-range(2,2)*Desc):(min(eledist)*DDensity):(range(2,2)*Desc));
z=((-range(2,3)*Desc):(min(eledist)*DDensity):(range(2,3)*Desc));
sizez=size(z,2);
vic_h=zeros(size(x,1),size(x,2),sizez);

figure
hold on
h=pcolor(x,y,vic_h(:,:,1)); 
h2=plot(node(:,1),node(:,2),'.k','MarkerSize',5);

for i=1:sizez -1

    node_current=node(((node(:,3)>=z(i))&(node(:,3)<z(i+1))),1:2);%映射到X-Y平面
    eledist_current=eledist(((node(:,3)>=z(i))&(node(:,3)<z(i+1))));
    %%        
    if isempty(node_current)
        vic_h(:,:,i)=-1;%-1作为标志
        if i==sizez-1
            vic_h(:,:,i)=-1;
        end    
    else
        for j=1:size(node_current,1)
            %计算影响系数矩阵
            temp=sqrt(  (x-node_current(j,1)).^2+...
                        (y-node_current(j,2)).^2);
            dis_cell=temp;
            temp(temp<=eledist_current(j))=eledist_current(j);    
            eff=Alpha.^(-(temp-eledist_current(j)));

            %计算尺寸函数
            si=eledist_current(j)-1+Beta.^temp;

            vic_h2=eff.*si;    
            vic_h(:,:,i)=vic_h(:,:,i)+vic_h2;

        end
%         vic_h(:,:,i)=vic_h(:,:,i)/(max(eledist)*j);
        vic_h(:,:,i)=vic_h(:,:,i)/(5*j);



        delete(h);
        delete(h2);
        h=pcolor(x,y,vic_h(:,:,i));
        h2=plot(node_current(:,1),node_current(:,2),'.k','MarkerSize',5);
        shading interp
        colormap(jet) 
        colorbar;

        pause(0.8);

            
    end
    fprintf('%d\n',i);
end
%%
%slice显示
% xslice = [-2.2,-1,0,1,2.2]; yslice =  [-2.2,-1,0,1,2.2]; zslice =  [-2.2,-1,0,1,2.2];
xslice = 0; yslice = [-1,1]; zslice = 0;
X=repmat(x,[1,1,sizez]);
Y=repmat(y,[1,1,sizez]);
Z=reshape(z,1,1,sizez);
Z=repmat(Z,[size(x,1),size(x,2),1]);
slice(X,Y,Z,vic_h,xslice*Amplify,yslice*Amplify,zslice*Amplify);
hold on;
plot3(node(:,1),node(:,2),node(:,3),'.k','MarkerSize',5);
colorbar;
colormap jet
axis equal;
xlabel ( '<--- X --->' );
ylabel ( '<--- Y --->' );
zlabel ( '<--- Z --->' );
title ( 'ply_filename' );

vic_h2=sum(vic_h2,3);
vic_h2=vic_h2/max(max(vic_h2));


% %显示
% vic_h2=vic_h2/max(max(vic_h2));
pcolor(x,y,vic_h(:,:,78));
shading interp
colormap(jet) 
colorbar;
hold on
plot(p(:,1),p(:,2),'.k','MarkerSize',5);

return
end


%%
%获取每个节点的实际尺寸
function eledist=elementdis(node,ele)
%1-2
d1=sqrt((node(ele(:,1),1)-node(ele(:,2),1)).^2+...
        (node(ele(:,1),2)-node(ele(:,2),2)).^2+...
        (node(ele(:,1),3)-node(ele(:,2),3)).^2);
%1-3
d2=sqrt((node(ele(:,1),1)-node(ele(:,3),1)).^2+...
        (node(ele(:,1),2)-node(ele(:,3),2)).^2+...
        (node(ele(:,1),3)-node(ele(:,3),3)).^2);
%2-3
d3=sqrt((node(ele(:,2),1)-node(ele(:,3),1)).^2+...
        (node(ele(:,2),2)-node(ele(:,3),2)).^2+...
        (node(ele(:,2),3)-node(ele(:,3),3)).^2);    

eleindex=ele(:,[1,2,1,3,2,3]);
k=sparse(eleindex,ones(size(eleindex)),ones(size(eleindex)));
k=full(k);
eledist=sparse(eleindex,ones(size(eleindex)),[d1,d1,d2,d2,d3,d3]);
eledist=full(eledist);
eledist=eledist./k;

return
end
