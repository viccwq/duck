%function MyTest(x,y)
function showTri3d(node,ele)
%%
%
%%
h=figure(2);
set(h,'rend','opengl');
hold on;
axis equal;
grid on;
trimesh(ele,node(:,1),node(:,2),node(:,3),'edgecolor','k','FaceAlpha',1,'FaceColor','w');
% plot3(node(:,1),node(:,2),node(:,3),'.r','MarkerSize',15); 
xlabel ( '<--- X --->' );
ylabel ( '<--- Y --->' );
zlabel ( '<--- Z --->' );
% title ( 'filename' );
view(3)
axis equal vis3d off;
rotate3d on
cameramenu;

return
end
