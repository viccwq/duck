%显示vol的slice 和 isosurface
function showDuck()
load duck
slice(permute(duck_vol,[2,1,3]),[60,110],[30],[45]);
% slice(permute(rat_vol,[2,1,3]),[50,150],[],[100]);
shading flat;
hold on
axis equal;
colorbar;
% h=trimesh(duck_element(:,1:3),duck_node(:,1),duck_node(:,2),duck_node(:,3),'facecolor','c','edgecolor','k','FaceAlpha',0.8);
xlabel('<--X--列>');
ylabel('<--Y--行>');
zlabel('<--Z-->');
end