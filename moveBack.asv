function [ output_point ] = moveBack(depth,ori_point,vec,out_nodes,norms)
%MOVEBACK Summary of this function goes here
%   Detailed explanation goes here
%function [ out_point ] = moveBack(depth,orin_point)
%将超出界限的节点按照其等值面的法向量，缩回到剖分区域
%out_point         m*3 返回后节点的坐标
%depth             三维深度矩阵
%orin_point        m*3 需要调整的点
%vec               ori_point将沿着该向量便平移
%%
%排除区域内的点
end_point=ori_point+vec;
output_point=zeros(size(ori_point));
[ outside,sizeVal ] = pointSize1_N(depth,end_point,size(depth));
outside=(sizeVal<-0.06);%-0.06与movePoint2.cpp中判定临界点的0.05一致
if any(outside)
    output_point(~outside,:)=end_point(~outside,:);%没有超出区域的点,先存储
    ori_node=ori_point(outside,:);%待处理的点的起点
    end_node=end_point(outside,:);%待处理的点的终点
    vec=vec(outside,:);           %待处理的点的增量
    sizeVal=sizeVal(outside);     %待处理的点的尺寸
    min(sizeVal)
    if ~isempty(find(sizeVal<-5))
        vec(sizeVal<-5,:)=cutTo(ori_node(sizeVal<-5,:),end_node(sizeVal<-5,:),vec(sizeVal<-5,:),depth,-5);
        end_node=ori_node+vec;
    end
    %转换成单位向量
    vec=vec./repmat(sqrt(sum(vec.^2,2)),1,3);
    
    %%
    %显示法向量
    % figure(2)
    % h1=plot3(end_node(:,1),end_node(:,2),end_node(:,3),'.r','MarkerSize',10); 
    %%
    %逐个缩回溢出点
    % fprintf('There are %d points to be moved\n',size(end_node,1));
    % node0=fix(end_node);
    % [tf, index] = ismember(node0,out_nodes,'rows');
    % norm0=norms(index,:);   
    % h2=line( [end_node(:,1)+3*norm0(:,1),end_node(:,1)]',...
    %         [end_node(:,2)+3*norm0(:,2),end_node(:,2)]',...
    %         [end_node(:,3)+3*norm0(:,3),end_node(:,3)]','Color','k');


    % for i=1:size(end_node,1)  
    %     node=ceil(end_node(i,:));
    %     if node(1)~=end_node(i,1)
    %         node(1)=node(1)-1;
    %     end
    %     if node(2)~=end_node(i,2)
    %         node(2)=node(2)-1;
    %     end
    %     if node(3)~=end_node(i,3)
    %         node(3)=node(3)-1;
    %     end
    %     node0=[ 0,0,0;
    %             0,1,0;
    %             1,0,0;
    %             1,1,0;
    %             0,0,1;
    %             0,1,1;
    %             1,0,1;
    %             1,1,1];
    %     node0=node0+repmat(node,8,1);
    %     node0=node;%直接查找该点
    %     [tf, index] = ismember(node0,out_nodes,'rows');
    %     
    %     index(index==0)=[];
    %     if isempty(index)
    %         pause;
    %     end
    %     norm0=sum(norms(index,:),1)/size(index,1);    
    %     end_node(i,:)=movePoint(depth,ori_node(i,:),end_node(i,:),norm0(i,:),vec(i,:));
    %     if rem(i,3)==0
    %         disp(i);
    %     end
    % end
    end_node=movePoint2(depth,size(depth),end_node,vec,out_nodes,norms);
%     delete(h1);
%     delete(h2);
    output_point(outside,:)=end_node;   
else
    output_point=end_point;
end

% p2 = patch(isosurface(xi,yi,zi,depth,0),...
%     'FaceColor','c','EdgeColor','none');
% isonormals(xi,yi,zi,depth,p2)
% view(3); axis tight equal
% camlight;  camlight(-80,-10); lighting phong; 
% title('Data Normals')


return
end

