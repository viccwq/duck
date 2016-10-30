function [node_new] = insertPoint(node,ele)
%INSERTPOINT Summary of this function goes here
%   Detailed explanation goes here
%在模型表面进行点的加密，添加三角形的重心和边的中点

%找出所有的边
pair=zeros(0,2);
localpairs=nchoosek(1:3,2); %排列组合1-2 1-3 1-4 2-3 ...
for ii=1:size(localpairs,1)     %获取拓扑结构的边
    pair=[pair;ele(:,localpairs(ii,:))];
end
pair=unique(sort(pair,2),'rows');
%中点
node_new=(node(pair(:,1),:)+node(pair(:,2),:))/2;
node_new=[node;node_new];
return
end

