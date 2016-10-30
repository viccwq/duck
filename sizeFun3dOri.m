function [ node_new,range,farthest ] = sizeFun3dOri(node)
%%
%SIZEFUN Summary of this function goes here
%最小边界点调整到原点
%[ node_new,range,farthest ] = sizeFun3d(node)
%输入：node        输入节点 n*3
%输出：node_new    调整后的节点坐标 n*3
%输出：range       调整后node在X Y Z方向上的最大值 1*3
%输出：farthest    离原点最远点的距离
%%
%获取点x,y,z 方向的范围,将图形的重心平移到原点
range=(minmax(node'))';
node_new=node-repmat(range(1,:),size(node,1),1);
range=range-repmat(range(1,:),2,1);
range(1,:)=[];
%离原点最远点的距离
temp=sqrt(sum(node_new(:,1).^2+node_new(:,2).^2+node_new(:,3).^2,2));
farthest=max(temp);
return
end

