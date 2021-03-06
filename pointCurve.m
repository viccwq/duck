function [ outPoint,outSize ] = pointCurve(node,ele,vol)
%POINTCURVE Summary of this function goes here
%计算节点的曲率
vol2=smooth3(vol,'box',7);
vol2=vol2.*(vol.^2);
vol2=smooth3(vol2,'box',5);


% % vec1=node(ele(:,2),:)-node(ele(:,1),:);
% % vec2=node(ele(:,3),:)-node(ele(:,1),:);
% % vecNormal=cross(vec1,vec2,2);
% % vecNormal=vecNormal./repmat(sqrt(sum(vecNormal.^2,2)),1,3);
center=(node(ele(:,1),:)+node(ele(:,2),:)+node(ele(:,3),:))/3;
% % endPoint=1*vecNormal+center;
% % [ temp,sizeVal ] = pointSize1_N(vol2,endPoint,size(vol));
% % if ~isempty(find(sizeVal==0))
% %     error('向量错误！');    
% % end
% % vecNormal(sizeVal>0,:)=vecNormal(sizeVal>0,:)*(-1);
% % endPoint=2*vecNormal+center;

[xi,yi,zi]=meshgrid(1:size(vol,1),1:size(vol,2),1:size(vol,3));
norms = isonormals(xi,yi,zi,permute(vol2,[2,1,3]),center);
norms=norms./repmat(sqrt(sum(norms.^2,2)),1,3);
endPoint=2*norms+center;
norms=repmat(norms,3,1);
% % isosurface(permute(vol2,[2,1,3]),0)
% % hold on
% % axis equal
% % showTri3d(node,ele);
% % hold on;
% % h=line([center(:,1),endPoint(:,1)]',...
% %        [center(:,2),endPoint(:,2)]',...
% %        [center(:,3),endPoint(:,3)]','Color','c','LineWidth',2);

pair=zeros(0,2);
localpairs=nchoosek(1:3,2); %排列组合1-2 1-3 1-4 2-3 ...
for ii=1:size(localpairs,1)     %获取拓扑结构的边
    pair=[pair;ele(:,localpairs(ii,:))];
end
pair=sort(pair,2);
[pair2,index1,temp]=unique(pair,'rows');
pair(index1,:)=0;
[tf,index2] = ismember(pair2, pair, 'rows');
tri=repmat([1:size(ele,1)]',3,1);
index1=tri(index1);%边邻接的一个三角形
index2=tri(index2);%边邻接的另一个三角形
pointSize=dot(norms(index1,:),norms(index2,:),2);
% % bars=pair2(find(pointSize<0.97),:);
% % h2=line([node(bars(:,1),1),node(bars(:,2),1)]',...
% %        [node(bars(:,1),2),node(bars(:,2),2)]',...
% %        [node(bars(:,1),3),node(bars(:,2),3)]','Color','r','LineWidth',2);
%% 边中点的size
outPoint=(node(pair2(:,1),:)+node(pair2(:,2),:))/2;
outSize=pointSize;
%% 表面点的size
pair2=pair2(:);
pointSize2=repmat(pointSize,2,1);
loc=full(sparse(pair2,ones(size(pair2)),pointSize2,size(node,1),1));
kk=full(sparse(pair2,ones(size(pair2)),ones(size(pair2)),size(node,1),1));
skip=kk==0;
outPoint2=node(~skip,:);
outSize2=loc(~skip)./kk(~skip);

outPoint=[outPoint;outPoint2];
outSize=[outSize;outSize2];
outSize=mapminmax(outSize',0,1)'*20;
outSize=2.2.^outSize;
% n=histc(outSize,0:max(outSize));
% bar(0:max(outSize),n,1,'c');
outSize=1-mapminmax(outSize',0,1)';
end

