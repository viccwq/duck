%function MyTest(x,y)
function showPoint3d(p)
%
%
%%


%%
%dol1
% p=p(p(:,2)>0,:);
%sphere
% p=p(p(:,1)>0&p(:,2)>0&p(:,3)>0,:);
h=figure(2);
set(h,'Renderer','OpenGL');
hold on;

plot3(p(:,1),p(:,2),p(:,3),'.r','MarkerSize',10); 
% tetramesh(t(num,:),p,'edgecolor','k','FaceAlpha',0.8,'FaceColor','g')
xlabel ( '<--- X --->' );
ylabel ( '<--- Y --->' );
zlabel ( '<--- Z --->' );
title ( 'point_filename' );

view(3)
axis equal vis3d off;
rotate3d on
cameramenu;

return
end
