function [ p,t ] = mesh3d( node,element )
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
[p,t]=distmeshnd_vic(fd,@huniform,2.0*mean(eledist),[-1,-1,-1;1,1,1]*radius,node);

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

if dim==1
    p=(box(1):h:box(2))';
elseif dim==3
    cbox=cell(1,dim);
    for ii=1:dim
        cbox{ii}=box(1,ii):h:box(2,ii);
    end
    pp=cell(1,dim);%笛卡尔网格
    [pp{:}]=ndgrid(cbox{:});
    p=zeros(prod(size(pp{1})),dim);
    for ii=1:dim
        p(:,ii)=pp{ii}(:);
    end
    
    %加密
%     load p2;    
%     p3=p2; p3(:,3)=p3(:,3)+0.2;
%     p4=p2; p4(:,3)=p4(:,2)+0.3;
%     p5=p2; p5(:,3)=p5(:,1)+0.2;

%cow_关闭加密
%     p3=p; p3(:,3)=p3(:,3)+0.4*5;
%     p4=p; p4(:,3)=p4(:,2)+0.3*5;
%     p5=p; p5(:,3)=p5(:,1)+0.7*5;
%     p=[p;p3;p4;p5];
%    

end

% 2. Remove points outside the region, apply the rejection method
p=p(feval(fdist,p,varargin{:})<geps,:);
%排除重复的点
p=p(~(uniquePoint(fix,p,6e-7)==1),:);
% p=p(~(uniquePoint(fix,p,2e-4)==1),:);

r0=feval(fh,p);
p=[fix; p(rand(size(p,1),1)<min(r0)^dim./r0.^dim,:)];   %概率选点

%*************************************************************************
 %构建八叉树
str='sphere';
str=[str,'_tree'];
 try
	load(str);    % Try block
    eval(['tree=',str,';']);
 catch ME
	[ fixInfor,tree ] = buildOctreeVic( box(:)',[fix,ones(size(fix,1),1)*1] ); 
    eval([str,'=tree;']);
    save(str,str);
 end
depth=pointDepth(tree,p);
p=p(depth~=0,:);
depth=depth(depth~=0);
%*************************************************************************

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
                simpplot(p,t,'p(:,1)>0');
%                 simpplot(p,t,'p(:,3)>0.4');
                title(['Retriangulation #',int2str(count)])
                drawnow
            end
        else
        	disp(sprintf('Retriangulation #%d',count))
        end
        count=count+1;
    end

  % 6. Move mesh points based on edge lengths L and forces F
  
  depth1=(depth(pair(:,1))+depth(pair(:,2)))/2;
  depth1=max(depth1)-depth1+1;
  depth1=depth1/max(depth1);
  
  bars=p(pair(:,1),:)-p(pair(:,2),:);
  L=sqrt(sum(bars.^2,2));                           %求边的长度
%   L0=feval(fh,(p(pair(:,1),:)+p(pair(:,2),:))/2);
%   L0=L0*L0mult*(sum(L.^dim)/sum(L0.^dim))^(1/dim);  %相对尺寸转化为本地尺寸 mean(L)<L0
%   F=max(L0-L,0);
%   F=max(L0-L,0).*(depth1+0.5);

    L0=depth1*L0mult*(sum(L.^dim)/sum(depth1.^dim))^(1/dim);
    F=L0-L;
    
    
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
