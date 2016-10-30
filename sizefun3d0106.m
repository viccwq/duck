function [ x,y,vic_h2 ] = sizefun3d( )
%SIZEFUN Summary of this function goes here
%   Detailed explanation goes here
%%
load cow_node;
load cow_element;
node=cow_node;
element=cow_element;
clear cow_node cow_element;

Desc=1.2;     %�ѿ��������С Descartes
DDensity=1/40;%�ѿ��������ܶȣ�Ϊͼ�γߴ�İٷֱ�
Alpha=20;      %Ӱ��ϵ������ ����
Beta=5;        %����ߴ纯�� ����
Amplify=10;    %����Ŵ�

node=node*Amplify;   %��������
eledist=elementdis(node,element);
% eledist=eledist/max(eledist);
%%
%��ȡ��ķ�Χ,��ͼ�ε�����ƽ�Ƶ�ԭ��
range=(minmax(node'))';
node=node-repmat((range(1,:)+range(2,:))/2,size(node,1),1);
range=range-repmat((range(1,:)+range(2,:))/2,2,1);
maxrange=max(range(2,:));
%�ѿ�������
[x,y]=meshgrid((-maxrange*Desc):(maxrange*DDensity):(maxrange*Desc));
z=((-maxrange*Desc):(maxrange*DDensity):(maxrange*Desc));
sizedesc=size(x,1); %size(x,1)==size(y,1)==size(z,1)
vic_h=zeros(sizedesc,sizedesc,sizedesc);
for i=1:sizedesc -1

    node_current=node(((node(:,3)>=z(i))&(node(:,3)<z(i+1))),1:2);%ӳ�䵽X-Yƽ��
    eledist_current=eledist(((node(:,3)>=z(i))&(node(:,3)<z(i+1))));
    %%
    if isempty(node_current)
        vic_h(:,:,i)=1/2;
        if i==sizedesc-1
            vic_h(:,:,i)=1/2;
        end    
    else
        for j=1:size(node_current,1)
            %����Ӱ��ϵ������
            temp=sqrt(  (x-node_current(j,1)).^2+...
                        (y-node_current(j,2)).^2);
            dis_cell=temp;
            temp(temp<=eledist_current(j))=eledist_current(j);    
            eff=Alpha.^(-(temp-eledist_current(j)));

            %����ߴ纯��
            si=eledist_current(j)-1+Beta.^temp;

            vic_h2=eff.*si;    
            vic_h(:,:,i)=vic_h(:,:,i)+vic_h2;
        end
        vic_h(:,:,i)=vic_h(:,:,i)/(max(eledist)*j);
    end
    fprintf('%d\n',i);
end
%%
%slice��ʾ
% xslice = [-2.2,-1,0,1,2.2]; yslice =  [-2.2,-1,0,1,2.2]; zslice =  [-2.2,-1,0,1,2.2];
xslice = 0; yslice = [-1,1]; zslice = 0;
X=repmat(x,[1,1,sizedesc]);
Y=repmat(y,[1,1,sizedesc]);
Z=reshape(z,1,1,sizedesc);
Z=repmat(Z,[sizedesc,sizedesc,1]);
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


% %��ʾ
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
%��ȡÿ���ڵ��ʵ�ʳߴ�
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
