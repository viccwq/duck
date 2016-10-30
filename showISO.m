close all
% [x y z v] = flow;
% [fc vt] = isosurface(x, y, z, v, -3);
load dol1_node;
load dol1_element;
vt=dol1_node;
fc=dol1_element;
clear dol1_node dol1_element;

p=patch('Faces',fc,'Vertices',vt);       
isonormals(x,y,z,v, p);       
set(p, 'FaceColor', 'w', 'EdgeColor', 'k');       
daspect([1 1 1]) ;      
view(3) ;     
camlight; 
% lighting phong;
% lighting gouraud£»
lighting flat
cameramenu;
figure;
plot3(vt(:,1),vt(:,2),vt(:,3),'.');