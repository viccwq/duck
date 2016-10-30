function showTree( tree )
%SHOWTREE Summary of this function goes here
%���ǰ˲���
%%
close all;
 tree=tree(((tree(:,1)>0)&(tree(:,2)>0)&(tree(:,3)>0)&(tree(:,4)>0)&(tree(:,5)>30)&(tree(:,6)>30)&(tree(:,5)<50)&(tree(:,6)<50)),:);
% tree=tree(((tree(:,1)>-10)&(tree(:,2)>-10)&(tree(:,3)>-10)&(tree(:,4)>-10)&(tree(:,5)>-0.3)&(tree(:,6)>-0.3)&(tree(:,5)<0.3)&(tree(:,6)<0.3)),:);
%tree=tree(((tree(:,1)>100)&(tree(:,2)>100)&(tree(:,3)>-10)&(tree(:,4)>-10)&(tree(:,5)>-10)&(tree(:,6)>-10)),:);
% tree=tree(((tree(:,1)>-10000)&(tree(:,2)>-10000)&(tree(:,3)>-000)&(tree(:,4)>-000)&(tree(:,5)>-1000)&(tree(:,6)>-1000)),:);
[m,n]=size(tree);
center=[sum(tree(:,1:2),2),sum(tree(:,3:4),2),sum(tree(:,5:6),2)]/2;




h=figure;
set(h,'Renderer','OpenGL');
hold on;

center1=center(tree(:,7)==7,:);
plot3(center1(:,1),center1(:,2),center1(:,3),'.r','MarkerSize',10); 
center1=center(tree(:,7)==6,:);
plot3(center1(:,1),center1(:,2),center1(:,3),'.g','MarkerSize',10); 
center1=center(tree(:,7)==5,:);
plot3(center1(:,1),center1(:,2),center1(:,3),'.b','MarkerSize',10); 
center1=center(tree(:,7)==4,:);
plot3(center1(:,1),center1(:,2),center1(:,3),'.k','MarkerSize',10); 
center1=center(tree(:,7)==3,:);
plot3(center1(:,1),center1(:,2),center1(:,3),'.m','MarkerSize',10); 
center1=center(tree(:,7)==2,:);
plot3(center1(:,1),center1(:,2),center1(:,3),'.r','MarkerSize',10); 
center1=center(tree(:,7)==1,:);
plot3(center1(:,1),center1(:,2),center1(:,3),'.c','MarkerSize',10); 
center1=center(tree(:,7)==0,:);
plot3(center1(:,1),center1(:,2),center1(:,3),'.y','MarkerSize',10); 
% text(center(:,1),center(:,2),center(:,3),[repmat('  ',m,1),num2str(tree(:,7))]); 
% tetramesh(t(num,:),p,'edgecolor','k','FaceAlpha',0.8,'FaceColor','g')
xlabel ( '<--- X --->' );
ylabel ( '<--- Y --->' );
zlabel ( '<--- Z --->' );
title ( 'Octree' );
view(3)
axis equal vis3d off;
rotate3d on

% cameramenu;

return
end

