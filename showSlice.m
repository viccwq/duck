function showSlice( tree )
%SHOWSLICE Summary of this function goes here
%   Detailed explanation goes here
close all;
a=minmax(tree(:,1:6)')';
b=max(a(2,:))*0.98;
[x,y,z]=meshgrid(-b:b*0.04:b);
str='sphere';
str=[str,'_depth'];

% try
%     load(str);
%     eval(['depth=',str,';']);
% catch ME
    depth= pointDepth( tree,[x(:),y(:),z(:)] );
    depth=reshape(depth,size(x));
%     eval([str,'=depth;']);
%     save(str,str);
% end

[xi,yi,zi]=meshgrid(-b:b*0.015:b);
vi = interp3(x,y,z,depth,xi,yi,zi);
% xslice = [-2.2,-1,0,1,2.2]; yslice =  [-2.2,-1,0,1,2.2]; zslice =  [-2.2,-1,0,1,2.2];
% xslice = [0.5,0]; yslice = [-1,1]; zslice = [0,0.4];
xslice = []; yslice = [0]; zslice = [0];
%�˲�
b= [1 1 1 1 1;
	1 1 1 1 1;
    1 1 1 1 1;
    1 1 1 1 1;
    1 1 1 1 1];
b1(:,:,1)=b;
b1(:,:,2)=b;
b1(:,:,3)=b; 
b1(:,:,4)=b;
b1(:,:,5)=b;
b1(2:4,2:4,2:4)=4;
b1(3,3,3)=4;
vi_filter=convn(vi,b1)/sum(sum(sum(b1,3),2));
vi_filter(1:2,:,:)=[];
vi_filter(:,1:2,:)=[];
vi_filter(:,:,1:2)=[];
vi_filter(end-1:end,:,:)=[];
vi_filter(:,end-1:end,:)=[];
vi_filter(:,:,end-1:end)=[];

slice(xi,yi,zi,vi_filter,xslice,yslice,zslice);
shading flat 
colorbar;
colormap jet;
axis equal;
xlabel ( '<--- X --->' );
ylabel ( '<--- Y --->' );
zlabel ( '<--- Z --->' );
title ( 'ply_filename' );
rotate3d on;
cameramenu;
end

