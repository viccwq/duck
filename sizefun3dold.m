function [ x,y,vic_h2 ] = sizefun3dold( )
%SIZEFUN Summary of this function goes here
%   Detailed explanation goes here
%%
Desc=1.1;     %笛卡尔网格大小 Descartes
DDensity=1/30;%笛卡尔网格密度，为图形尺寸的百分比
Alpha=5;      %影响系数矩阵 底数
Beta=8;        %%计算尺寸函数 底数
Amplify=10;    %坐标放大
Limit=1e-7;    %为0的阈值
load cow_node;
load cow_element;
cow_node=cow_node*Amplify;
eledist=elementdis(cow_node,cow_element);
% eledist=eledist/max(eledist);
%%
%获取点的范围,将图形的重心平移到原点
range=(minmax(cow_node'))';
cow_node=cow_node-repmat((range(1,:)+range(2,:))/2,size(cow_node,1),1);
range=range-repmat((range(1,:)+range(2,:))/2,2,1);
maxrange=max(range(2,:));
%笛卡尔网格
[x,y,z]=meshgrid((-maxrange*Desc):(maxrange*DDensity):(maxrange*Desc));
sizedesc=size(x,1);
vic_h=zeros(sizedesc,sizedesc,sizedesc);
for i=1:size(cow_node,1)
%     %计算影响系数矩阵
%     temp=sqrt(  (x-ones(sizedesc,sizedesc,sizedesc)*cow_node(i,1)).^2+...
%                 (y-ones(sizedesc,sizedesc,sizedesc)*cow_node(i,2)).^2+...
%                 (z-ones(sizedesc,sizedesc,sizedesc)*cow_node(i,3)).^2);
%     dis_cell{i}=temp;
%     temp(temp<=eledist(i))=eledist(i);    
%     eff{i}=Alpha.^(-(temp-eledist(i)));
%     
%     %计算尺寸函数
%     si{i}=dis_cell{i}-1+Beta.^temp;
%     
%     vic_h2{i}=eff{i}.*si{i};    
%     vic_h=vic_h+vic_h2{i};
%     fprintf('%d\n',i);

    %计算影响系数矩阵
    temp=sqrt(  (x-ones(sizedesc,sizedesc,sizedesc)*cow_node(i,1)).^2+...
                (y-ones(sizedesc,sizedesc,sizedesc)*cow_node(i,2)).^2+...
                (z-ones(sizedesc,sizedesc,sizedesc)*cow_node(i,3)).^2);
    dis_cell=temp;
    temp(temp<=eledist(i))=eledist(i);    
    eff=Alpha.^(-(temp-eledist(i)));
    eff(eff<Limit)=0;
    %计算尺寸函数
%     si=eledist(i)-1+Beta.^temp;
    si=(1-Beta.^(-temp))*2+eledist(i);

    vic_h2=eff.*si;    
    vic_h=vic_h+vic_h2;
    fprintf('%d\n',i);



end
%%
% vic_h2=sum(vic_h2,3);
% vic_h2=vic_h2/max(max(vic_h2));
vic_h=vic_h/max(max(max(vic_h,[],3)));


%slice显示
% xslice = [-2.2,-1,0,1,2.2]; yslice =  [-2.2,-1,0,1,2.2]; zslice =  [-2.2,-1,0,1,2.2];
% xslice = 0; yslice = [-1,1]; zslice = 0;
% slice(x,y,z,vic_h,xslice*Amplify,yslice*Amplify,zslice*Amplify);
xslice = cow_node(i,1); yslice = cow_node(i,2); zslice = cow_node(i,3);
slice(x,y,z,vic_h,xslice ,yslice ,zslice );
hold on;
plot3(cow_node(:,1),cow_node(:,2),cow_node(:,3),'.k','MarkerSize',5);
colorbar;
colormap jet;
axis equal;
xlabel ( '<--- X --->' );
ylabel ( '<--- Y --->' );
zlabel ( '<--- Z --->' );
title ( 'ply_filename' );
rotate3d on;



% %显示
% vic_h2=vic_h2/max(max(vic_h2));
pcolor(x,y,vic_h2);
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
