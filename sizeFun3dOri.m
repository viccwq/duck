function [ node_new,range,farthest ] = sizeFun3dOri(node)
%%
%SIZEFUN Summary of this function goes here
%��С�߽�������ԭ��
%[ node_new,range,farthest ] = sizeFun3d(node)
%���룺node        ����ڵ� n*3
%�����node_new    ������Ľڵ����� n*3
%�����range       ������node��X Y Z�����ϵ����ֵ 1*3
%�����farthest    ��ԭ����Զ��ľ���
%%
%��ȡ��x,y,z ����ķ�Χ,��ͼ�ε�����ƽ�Ƶ�ԭ��
range=(minmax(node'))';
node_new=node-repmat(range(1,:),size(node,1),1);
range=range-repmat(range(1,:),2,1);
range(1,:)=[];
%��ԭ����Զ��ľ���
temp=sqrt(sum(node_new(:,1).^2+node_new(:,2).^2+node_new(:,3).^2,2));
farthest=max(temp);
return
end

