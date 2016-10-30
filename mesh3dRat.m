function [ p,t ] = mesh3dRat( node,element,vol)
%%
%20130302 ������
%20130307 ���Ժ���
%20130311 smooth3
%20130313 rat
%%
%MESH3D Summary of this function goes here
%������ά����ģ�ͣ�����delaunay�����ʷ�
%[ output  ] = mesh3d( node,element )
%node��      ����ڵ����� n*3
%element:    ���������ε����� m*3


%%
%��������
% Desc=1.2;     %�ѿ��������С Descartes
% DDensity=2;%�ѿ��������ܶȣ�Ϊͼ����С�ߴ�İٷֱ�
% Alpha=20;      %Ӱ��ϵ������ ����
% Beta=12;        %����ߴ纯�� ����
Amplify=1;    %����Ŵ�
swith=false;
%%
%�������
if swith
node=node*Amplify;      %�Ŵ�����
[ node_new,range,farthest ] = sizeFun3d(node);
node=node_new;
clear node_new;
end
%%
rand('state',1); % Always the same results
set(gcf,'rend','opengl');
disp('Rat head');
radius=0;
fd=inline(strcat('sqrt(sum(p.^2,2))-',num2str(radius)),'p');
% [p,t]=distmeshnd_vic(fd,@huniform,0.2*radius,[-1,-1,-1;1,1,1]*radius,node);
%�ų�surface���ܵĵ�
eledist=nodeSize(node,element);
fix=node(~uniquePointSelf (node,mean(eledist)*0.35),:);

node=insertPoint(node,element);%����ڵ����
node=node(~uniquePointSelf (node,mean(eledist)*0.3),:);

vol=volConvert(vol);
box=[0,0,0;size(vol)+1];
[p,t]=distmeshnd_vic(fd,@huniform,mean(eledist),box,fix,node,vol);

% MyTest3d(p);
save pANDt3d p t;
% post(p,t)


return
end


%%
function post(p,t)

disp(sprintf('(press any key)'))
disp('����')

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
function vol2=volConvert(vol,varargin)
count=size(varargin,2);
if count==0
    %��vol��׼����ת��Ϊint8
    vol(vol==0)=vol(vol==0)+3;
    vol2=vol*(-1)+2;
else
    %vol��depth���
    depth=varargin{1};
    if size(vol)==size(depth)
        vol2=depth.*double(vol);
    else
        vol=permute(vol,[2,1,3]);
        if size(vol)==size(depth)
            vol2=depth.*double(vol);
        else
            error('vol��depth�Ĵ�С��һ��');
        end
    end
end

return
end
%%
function [p,t]=distmeshnd_vic(fdist,fh,h,box,fix,node,vol,varargin)
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
[xi,yi,zi]=meshgrid(1:box(2,1)-1,1:box(2,2)-1,1:box(2,3)-1);
% [x,y,z]=meshgrid(1:2:box(2,1)-1,1:2:box(2,2)-1,1:2:box(2,3)-1);
% clear xi yi zi
dim=size(box,2);

ptol=.001; 
ttol=.05; %0.1
L0mult=1+.4/2^(dim-1);  %�����ߴ�ķŴ�
deltat=.1;              %�仯ʱ��
geps=1e-1*h;            %���ޣ��ų������
deps=sqrt(eps)*h;

% 1. Create initial distribution in bounding box

swith=false;
if swith
%*************************************************************************
%�����˲���

%�������ߴ�
stro='rat2';
str=[stro,'_tree'];
 try
	load(str);    % Try block
    eval(['tree=',str,';']);
 catch ME
	[ fixInfor,tree ] = buildOctreeVic( box(:)',[node,ones(size(node,1),1)*7] ); 
    eval([str,'=tree;']);
    save(str,str);
 end
eval(['clear ',str,';']);
clear fixInfor
%*************************************************************************


%������С��Ԫ
% mintree=find(tree(:,7)==max(tree(:,7)));
% mintree=tree(mintree(1),2)-tree(mintree(1),1);


%�����������

str=[stro,'_depth'];
try
    load(str);
    eval(['depth=',str,';']);
catch ME
    depth= pointDepth( tree,[xi(:),yi(:),zi(:)] );
    depth=reshape(depth,size(xi));
    eval([str,'=depth;']);
    save(str,str);
end
eval(['clear ',str,';']);

clear tree

%����������ֵ
% depthMax=max(max(max(depth,[],3)));
% depth(depth==depthMax)=depth(depth==depthMax)-1;
% clear depthMax;
%%
%�˲�1
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

clear de

% depth=max(max(max(depth,[],3)))+min(min(min(depth,[],3)))-depth;
depth=max(max(max(depth,[],3)))-depth+1;
depth=smooth3(depth,'box',3);
depth=log(depth)/log(1.2);


%%
%�˲�2
depth=depth-min(min(min(depth,[],3)));
depth=depth.*permute(vol.^2,[2,1,3]);
depth=smooth3(depth,'box',3);
depth=smooth3(depth,'box',5);
depth=depth.*permute(vol.^2,[2,1,3]);
depth=smooth3(depth,'box',3);
depth=depth/max(max(max(depth,[],3)));
depth=depth*4;
depth=depth-min(min(min(depth,[],3)))+1;
depth=log(depth)/log(1.2);
depth=smooth3(depth,'box',5);
depth=depth-min(min(min(depth,[],3)));
depth=(1.2-0.05).^depth;
depth=depth-min(min(min(depth,[],3)))+0.05;
depth=depth.*permute(vol,[2,1,3]);
depth=smooth3(depth,'box',3);



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
end




% ��ʾ��Ƭ
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


%*************************************************************************
clear vol
load depth
depth=depth/max(max(max(abs(depth),[],3)))*10;
% ��ʾ��ֵ��
% figure(2)
% isosurface(depth,0);
% axis equal
% view(3)
% hold on

%%
p=[xi(:),yi(:),zi(:)];
[ temp,sizeVal ] = pointSize1_N(depth,p,size(depth));

%���������ķ�����
out_nodes=p(sizeVal<4&sizeVal>-5,:);
norms = isonormals(xi,yi,zi,depth,out_nodes);
%ת���ɵ�λ����
norms=norms./repmat(sqrt(sum(norms.^2,2)),1,3);
clear xi yi zi ;
% 2. Remove points outside the region, apply the rejection method
p=p(sizeVal>0.05,:);%0.2
sizeVal=sizeVal(sizeVal>0.05);%0.2
%ƽ������ɸѡ
temp=rand(size(p,1),1)<0.26;%0.2���� 0.25�е���̬
p=p(temp,:);
sizeVal=sizeVal(temp);
clear temp
r0=sizeVal/max(sizeVal);
p=p(rand(size(p,1),1)<(min(r0)+0.05)^(3/2)./r0.^(3/2),1:3);
[ temp ] = uniquePointSelf (p,1.5);
p=p(temp==0,:);

[ outData ] = uniquePointSelf (fix,5);
fix=fix(outData==0,:);
%�ų��ظ��ĵ�
p=p(~(uniquePoint(fix,p(:,1:3),6e-7)==1),:);
% p=[fix;p];


N=size(p,1);
count=0;
p0=inf;%�������ȡ���λ��
while 1
  % 3. Retriangulation by Delaunay
    figure(1)
    %���λ�Ƴ���һ��ֵ�����½������˽ṹ
    if max(sqrt(sum((p-p0).^2,2)))>ttol*h   
        p0=p;
        t=delaunayn(p); %��Ҫ�����ظ��ĵ㣡����������������
        pmid=zeros(size(t,1),dim);
        for ii=1:dim+1  %�����������
            pmid=pmid+p(t(:,ii),:)/(dim+1);
        end
        [ temp,sizeVal ] = pointSize1_N(depth,pmid,size(depth));
        t=t(sizeVal>=0.2,:); %�ų������������������
        % 4. Describe each edge by a unique pair of nodes
        pair=zeros(0,2);
        localpairs=nchoosek(1:dim+1,2); %�������1-2 1-3 1-4 2-3 ...
        for ii=1:size(localpairs,1)     %��ȡ���˽ṹ�ı�
            pair=[pair;t(:,localpairs(ii,:))];
        end
        pair=unique(sort(pair,2),'rows');
        % 5. Graphical output of the current mesh
        if dim==2
            trimesh(t,p(:,1),p(:,2),zeros(N,1))
            view(2),axis equal,axis off,drawnow
        
        elseif dim==3
            if mod(count,10)==0
                simpplot(p,t,'p(:,2)>50');
%                 simpplot(p,t,'p(:,3)>0.4');
                title(['Retriangulation #',int2str(count)])
                drawnow
            end
        else
        	disp(sprintf('Retriangulation #%d',count))
        end
        count=count+1;
    end
    
    if rem(count,6)==0
        showStaQ(p,t);
    end
    
    

  % 6. Move mesh points based on edge lengths L and forces F
  bars=p(pair(:,1),:)-p(pair(:,2),:);
  L=sqrt(sum(bars.^2,2));                           %��ߵĳ���
  
  [ temp,sizeVal ] = pointSize1_N(depth,(p(pair(:,1),:)+p(pair(:,2),:))/2,size(depth));
  depth1=sizeVal;
  depth1(sizeVal<1)=1;  
  depth1=depth1./depth1;

  L0=depth1*L0mult*(sum(L.^(dim-1))/sum(depth1.^(dim-1)))^(1/(dim-1));
%     F=L0-L;
  F=max(L0-L,0);    
  Fbar=[bars,-bars].*repmat(F./L,1,2*dim);          %bar����[x,y,x].*repmat(F./L,1,3);  [-x,-y,-z].*repmat(F./L,1,3)
  dp=full(sparse(pair(:,[ones(1,dim),2*ones(1,dim)]), ...%��ÿ��ĺ���
                 ones(size(pair,1),1)*[1:dim,1:dim], ...
                 Fbar,N,dim));
  dp(1:size(fix,1),:)=0;                            %�̶������Ϊ0
%   p=p+deltat*dp;

  % 7. Bring outside points back to the boundary
  %�������
%   p(size(fix,1)+1:end,:) = moveBack(depth,p(size(fix,1)+1:end,:),deltat*dp(size(fix,1)+1:end,:),out_nodes,norms);%out_nodes��norms������ѯ��ķ�����
  %�ų������
  p= moveBack(depth,p,deltat*dp,out_nodes,norms);%out_nodes��norms������ѯ��ķ�����
  
  
  if swith
  d=feval(fdist,p,varargin{:}); ix=d>0;
  gradd=zeros(sum(ix),dim);
  for ii=1:dim
    a=zeros(1,dim);
    a(ii)=deps;
    d1x=feval(fdist,p(ix,:)+ones(sum(ix),1)*a,varargin{:});
    gradd(:,ii)=(d1x-d(ix))/deps;
  end
  p(ix,:)=p(ix,:)-d(ix)*ones(1,dim).*gradd;
  end
  % 8. Termination criterion
%   maxdp=max(deltat*sqrt(sum(dp(d<-geps,:).^2,2)));
%   if maxdp<ptol*h, break; end
  if count==100, break; end
end
return
end
