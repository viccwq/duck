function [node_new] = insertPoint(node,ele)
%INSERTPOINT Summary of this function goes here
%   Detailed explanation goes here
%��ģ�ͱ�����е�ļ��ܣ���������ε����ĺͱߵ��е�

%�ҳ����еı�
pair=zeros(0,2);
localpairs=nchoosek(1:3,2); %�������1-2 1-3 1-4 2-3 ...
for ii=1:size(localpairs,1)     %��ȡ���˽ṹ�ı�
    pair=[pair;ele(:,localpairs(ii,:))];
end
pair=unique(sort(pair,2),'rows');
%�е�
node_new=(node(pair(:,1),:)+node(pair(:,2),:))/2;
node_new=[node;node_new];
return
end

