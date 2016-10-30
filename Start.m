clc;
close all;
[xi,yi,zi]=meshgrid(1:221,1:209,1:201);
[x,y,z]=meshgrid(1:2:221,1:2:209,1:2:201);
clear xi yi zi x y z
%%
%cow
% load cow_node;
% load cow_element;
% node=cow_node;
% element=cow_element;
% clear cow_node cow_element;
%%
%sphere
% load sphere_node;
% load sphere_element;
% node=sphere_node;
% element=sphere_element;
% clear sphere_node sphere_element;
%%
% dol1
% load dol1_node;
% load dol1_element;
% node=dol1_node;
% element=dol1_element;
% clear dol1_node dol1_element;
%%
%rat
% load rat;
% node=rat_node;
% element=rat_face;
% vol=rat_vol;
% clear rat_node rat_face rat_vol;
%%
%rat2
% load rat2;
% node=rat2_node;
% element=rat2_face;
% vol=rat2_vol;
% clear rat2_node rat2_face rat2_vol;
%%
%duck
load duck;
node=duck_node;
element=duck_element;
vol=duck_vol;
clear duck_node duck_element duck_vol

%%
% recoverVol( node,element);
mesh3dDuck(node,element,vol);
