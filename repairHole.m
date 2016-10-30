function [] = repairHole(node,surf,addNode)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
clc
node=[node;addNode];
h=trimesh(surf,node(:,1),node(:,2),node(:,3),'facecolor',[.001,1,1],'edgecolor','y','FaceAlpha',0.5);
hold on
pair=zeros(0,2);
localpairs=nchoosek(1:3,2); %排列组合1-2 1-3 1-4 2-3 ...
for ii=1:size(localpairs,1)     %获取拓扑结构的边
    pair=[pair;surf(:,localpairs(ii,:))];
end
pair=sort(pair,2);
[bars,ix,jx]=unique(pair,'rows');
vec=histc(jx,1:max(jx));
qx=find(vec==1);
unibar=pair(ix(qx),:);
line(   [node(unibar(:,1),1),node(unibar(:,2),1)]',...
        [node(unibar(:,1),2),node(unibar(:,2),2)]',...
        [node(unibar(:,1),3),node(unibar(:,2),3)]','Color','r','LineWidth',3);
hold on
% plot3(node(unibar(:,1),1),node(unibar(:,1),2),node(unibar(:,1),3),'r');
% scatter3(node(unibar(:,1),1),node(unibar(:,1),2),node(unibar(:,1),3),ones(size(node(unibar(:,1),3))),'filled');
unibar=unique(unibar(:));
scatter3(node(unibar,1),node(unibar,2),node(unibar,3),ones(size(node(unibar(:,1),3)))*5,'filled');
txtpars={'fontname','times','fontsize',14,'horizontala','left'};
text(node(unibar,1),node(unibar,2),node(unibar,3),[repmat('    ',size(unibar)),num2str(unibar)],txtpars{:});
if ~isempty(addNode)
    scatter3(addNode(:,1),addNode(:,2),addNode(:,3),ones(size(addNode,1),1)*15,'filled');
end
view(3)
axis equal vis3d;
rotate3d on
grid off
% cameramenu;

end

