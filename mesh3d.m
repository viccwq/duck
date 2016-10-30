function [ p,t ] = mesh3d( node,element )
%%
%20130302 调试球
%20130307 调试海豚
%20130311 smooth3
%%
%MESH3D Summary of this function goes here
%输入三维表面模型，进行delaunay三角剖分
%[ output  ] = mesh3d( node,element )
%node：      输入节点坐标 n*3
%element:    表面三角形的索引 m*3


%%
%参数设置
% Desc=1.2;     %笛卡尔网格大小 Descartes
% DDensity=2;%笛卡尔网格密度，为图形最小尺寸的百分比
% Alpha=20;      %影响系数矩阵 底数
% Beta=12;        %计算尺寸函数 底数
Amplify=1;    %坐标放大

%%
%坐标调整
node=node*Amplify;      %放大坐标
[ node_new,range,farthest ] = sizeFun3d(node);
node=node_new;
clear node_new;

%%
rand('state',1); % Always the same results
set(gcf,'rend','opengl');

disp('(9) 3-D Unit ball');
radius=farthest*1.2;
fd=inline(strcat('sqrt(sum(p.^2,2))-',num2str(radius)),'p');
% [p,t]=distmeshnd_vic(fd,@huniform,0.2*radius,[-1,-1,-1;1,1,1]*radius,node);
eledist=nodeSize(node,element);

node=insertPoint(node,element);%表面节点加密
[p,t]=distmeshnd_vic(fd,@huniform,mean(eledist),[-1,-1,-1;1,1,1]*radius,node);

% MyTest3d(p);
save pANDt3d p t;
% post(p,t)


return
end


%%
function post(p,t)

disp(sprintf('(press any key)'))
disp('结束')

end

%%
function d=fd10(p)

r=sqrt(p(:,1).^2+p(:,2).^2);
z=p(:,3);

d1=r-1;
d2=z-1;
d3=-z-1;
d4=sqrt(d1.^2+d2.^2);
d5=sqrt(d1.^2+d3.^2);
d=dintersect(dintersect(d1,d2),d3);
ix=d1>0 & d2>0;
d(ix)=d4(ix);
ix=d1>0 & d3>0;
d(ix)=d5(ix);

d=ddiff(d,dsphere(p,0,0,0,0.5));
return
end


%%
function h=fh10(p)

h1=4*sqrt(sum(p.^2,2))-1;
h=min(h1,2);
return
end


%%
function [p,t]=distmeshnd_vic(fdist,fh,h,box,fix,varargin)
%DISTMESHND N-D Mesh Generator using Distance Functions.
%   [P,T]=DISTMESHND(FDIST,FH,H,BOX,FIX,FDISTPARAMS)
%
%      P:           Node positions (NxNDIM)
%      T:           Triangle indices (NTx(NDIM+1))
%      FDIST:       Distance function
%      FH:          Edge length function
%      H:           Smallest edge length
%      BOX:         Bounding box [xmin,xmax;ymin,ymax; ...] (NDIMx2)
%      FIX:         Fixed node positions (NFIXxNDIM)
%      FDISTPARAMS: Additional parameters passed to FDIST
%
%   Example: Unit ball
%      dim=3;
%      d=inline('sqrt(sum(p.^2,2))-1','p');
%      [p,t]=distmeshnd(d,@huniform,0.2,[-ones(1,dim);ones(1,dim)],[]);
%
%   See also: DISTMESH2D, DELAUNAYN, TRIMESH, MESHDEMOND.

%   Copyright (C) 2004-2012 Per-Olof Persson. See COPYRIGHT.TXT for details.

dim=size(box,2);

ptol=.001; 
ttol=.1; 
L0mult=1+.4/2^(dim-1);  %期望尺寸的放大
deltat=.1;              %变化时间
geps=1e-1*h;            %下限，排除球外点
deps=sqrt(eps)*h;

% 1. Create initial distribution in bounding box


%*************************************************************************
 %构建八叉树
stro='dol1';
str=[stro,'_tree'];
 try
	load(str);    % Try block
    eval(['tree=',str,';']);
 catch ME
	[ fixInfor,tree ] = buildOctreeVic( box(:)',[fix,ones(size(fix,1),1)*7] ); 
    eval([str,'=tree;']);
    save(str,str);
 end

%*************************************************************************


%树的最小单元
mintree=find(tree(:,7)==max(tree(:,7)));
mintree=tree(mintree(1),2)-tree(mintree(1),1);

bound=box(2,1)*0.95;
ans=bound/h;

%计算各点的深度
%dol1
[x,y,z]=meshgrid(-bound:bound/(ceil(ans)+45):bound);

%sphere
% [x,y,z]=meshgrid(-bound:bound/(ceil(ans)*1):bound);

str=[stro,'_depth'];
try
    load(str);
    eval(['depth=',str,';']);
catch ME
    depth= pointDepth( tree,[x(:),y(:),z(:)] );
    depth=reshape(depth,size(x));
    eval([str,'=depth;']);
    save(str,str);
end
%降低最高深度值
% depthMax=max(max(max(depth,[],3)));
% depth(depth==depthMax)=depth(depth==depthMax)-1;
% clear depthMax;
%%
%滤波1
de=zeros(size(depth)+2);
depth=[depth(1,:,:);depth;depth(end,:,:)];
depth=[depth(:,1,:),depth,depth(:,end,:)];
de(:,:,1)=depth(:,:,1);
de(:,:,2:end-1)=depth;
de(:,:,end)=depth(:,:,end);

b(:,:,1)= [1 1 1;1 2 1;1 1 1];
b(:,:,2)= [1 2 1;2 2 2;1 2 1];
b(:,:,3)= [1 1 1;1 2 1;1 1 1];
% b(:,:,1)= [2 2 2;2 1 2;2 2 2];
% b(:,:,2)= [2 1 2;1 1 1;2 1 2];
% b(:,:,3)= [2 2 2;2 1 2;2 2 2];
% b(:,:,1)= [1 1 1;1 1 1;1 1 1];
% b(:,:,2)= [1 1 1;1 1 1;1 1 1];
% b(:,:,3)= [1 1 1;1 1 1;1 1 1];
depth=convn(de,b,'same')/sum(sum(sum(b,3),2));

depth(1,:,:)=[];
depth(:,1,:)=[];
depth(:,:,1)=[];
depth(end,:,:)=[];
depth(:,end,:)=[];
depth(:,:,end)=[];

% depth=max(max(max(depth,[],3)))+min(min(min(depth,[],3)))-depth;
depth=max(max(max(depth,[],3)))-depth+1;
depth=smooth3(depth,'box',3);
depth=log(depth)/log(1.2);



%%
%滤波2
%dol1
[xi,yi,zi]=meshgrid(-bound:bound/(ceil(ans)+60):bound);
%sphere
% [xi,yi,zi]=meshgrid(-bound:bound/(ceil(ans)*6):bound);
depth = interp3(x,y,z,depth,xi,yi,zi);

depth=smooth3(depth,'box',5);


% de=zeros(size(depth)+4);
% depth=[depth(1:2,:,:);depth;depth(end-1:end,:,:)];
% depth=[depth(:,1:2,:),depth,depth(:,end-1:end,:)];
% de(:,:,1:2)=depth(:,:,1:2);
% de(:,:,3:end-2)=depth;
% de(:,:,end-1:end)=depth(:,:,end-1:end);
% b= [1 1 1 1 1;
% 	1 1 1 1 1;
%     1 1 1 1 1;
%     1 1 1 1 1;
%     1 1 1 1 1];
% b=repmat(b,[1,1,5]);
% b(2:4,2:4,2:4)=0.5;
% b(3,3,3)=0.2;
% depth=convn(de,b,'same')/sum(sum(sum(b,3),2));
% depth(1:2,:,:)=[];
% depth(:,1:2,:)=[];
% depth(:,:,1:2)=[];
% depth(end-1:end,:,:)=[];
% depth(:,end-1:end,:)=[];
% depth(:,:,end-1:end)=[];



% depth=1.8.^depth;
depth=1.3.^depth;
depth=depth-min(min(min(depth,[],3)))*0.8;

% 显示切片
if 0
    figure(2);
    xslice = [0]; yslice = [0]; zslice = [0];
    slice(xi,yi,zi,depth,xslice,yslice,zslice);
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
scale=xi(1,:,1);
clear x y z de b str;
%*************************************************************************
%%
p=[xi(:),yi(:),zi(:),depth(:)];
%平均概率筛选
% p=p(rand(size(p,1),1)<0.7,:);


figure(1);
% 2. Remove points outside the region, apply the rejection method
p=p(feval(fdist,p(:,1:3),varargin{:})<geps,:);
%排除重复的点
p=p(~(uniquePoint(fix,p(:,1:3),6e-7)==1),:);
% p=p(~(uniquePoint(fix,p,2e-4)==1),:);

r0=p(:,4);

% r0=r0/max(r0);
%dol1
fix(286:end,:)=[];
p=[fix; p(rand(size(p,1),1)<(min(r0)+0.18)^3./r0.^3,1:3)];%概率选点
% p= p(rand(size(p,1),1)<(min(r0)+0.0)^3./r0.^3,1:3);
%sphere
% p=[fix; p(rand(size(p,1),1)<(min(r0)+0.0)^5./r0.^5,1:3)];   %概率选点



N=size(p,1);
count=0;
p0=inf;%正无穷，获取最大位移
while 1
  % 3. Retriangulation by Delaunay
    %最大位移超过一定值，重新建立拓扑结构
    if max(sqrt(sum((p-p0).^2,2)))>ttol*h   
        p0=p;
        t=delaunayn(p); %需要处理重复的点！！！！！！！！！
        pmid=zeros(size(t,1),dim);
        for ii=1:dim+1  %四面体的重心
            pmid=pmid+p(t(:,ii),:)/(dim+1);
        end
        t=t(feval(fdist,pmid,varargin{:})<-geps,:); %排除重心在球外的四面体
        % 4. Describe each edge by a unique pair of nodes
        pair=zeros(0,2);
        localpairs=nchoosek(1:dim+1,2); %排列组合1-2 1-3 1-4 2-3 ...
        for ii=1:size(localpairs,1)     %获取拓扑结构的边
            pair=[pair;t(:,localpairs(ii,:))];
        end
        pair=unique(sort(pair,2),'rows');
        % 5. Graphical output of the current mesh
        if dim==2
            trimesh(t,p(:,1),p(:,2),zeros(N,1))
            view(2),axis equal,axis off,drawnow
        
        elseif dim==3
            if mod(count,10)==0
                simpplot(p,t,'p(:,2)>0');
%                 simpplot(p,t,'p(:,3)>0.4');
                title(['Retriangulation #',int2str(count)])
                drawnow
            end
        else
        	disp(sprintf('Retriangulation #%d',count))
        end
        count=count+1;
    end
    
    if rem(count,36)==0
        showStaQ(p,t);
    end
    
    

  % 6. Move mesh points based on edge lengths L and forces F
  
  depth1=pointSize(scale',depth,(p(pair(:,1),:)+p(pair(:,2),:))/2);

%   depth1=depth1/max(depth1);
  depth1=depth1;
  
  bars=p(pair(:,1),:)-p(pair(:,2),:);
  L=sqrt(sum(bars.^2,2));                           %求边的长度
%   L0=feval(fh,(p(pair(:,1),:)+p(pair(:,2),:))/2);
%   L0=L0*L0mult*(sum(L.^dim)/sum(L0.^dim))^(1/dim);  %相对尺寸转化为本地尺寸 mean(L)<L0
%   F=max(L0-L,0);
%   F=max(L0-L,0).*(depth1+0.5);

  L0=depth1*L0mult*(sum(L.^dim)/sum(depth1.^dim))^(1/dim);
%     F=L0-L;
  F=max(L0-L,0);
    
    
  Fbar=[bars,-bars].*repmat(F./L,1,2*dim);          %bar向量[x,y,x].*repmat(F./L,1,3);  [-x,-y,-z].*repmat(F./L,1,3)
  dp=full(sparse(pair(:,[ones(1,dim),2*ones(1,dim)]), ...%求每点的合力
                 ones(size(pair,1),1)*[1:dim,1:dim], ...
                 Fbar,N,dim));
  dp(1:size(fix,1),:)=0;                            %固定点的力为0
  p=p+deltat*dp;

  % 7. Bring outside points back to the boundary
  d=feval(fdist,p,varargin{:}); ix=d>0;
  gradd=zeros(sum(ix),dim);
  for ii=1:dim
    a=zeros(1,dim);
    a(ii)=deps;
    d1x=feval(fdist,p(ix,:)+ones(sum(ix),1)*a,varargin{:});
    gradd(:,ii)=(d1x-d(ix))/deps;
  end
  p(ix,:)=p(ix,:)-d(ix)*ones(1,dim).*gradd;

  % 8. Termination criterion
  maxdp=max(deltat*sqrt(sum(dp(d<-geps,:).^2,2)));
  if maxdp<ptol*h, break; end
end
return
end
