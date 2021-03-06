function [ num,inside ] = inTetra( p,tetra,node )
%% 
%INTETRA
%判断点集node所在的四面体单元的序列号
%输入:p       输入点集
%输入:tetra   四面体（一行四列）,和p构成四面体网格
%输入:node    需要检测的点集 m*n
%输出：num    node中每个点所在的四面体的序列 m*1
%inside =1  在四面体内部 m*1
%       =0  在四面体的面、棱、顶点
%       =-1 在四面体外
%% 
%web http://steve.hollasch.net/cgindex/geometry/ptintet.html

%%
%dzero=1.45e-14;
dzero=1.45e-9;
pt=node;      %待检测的点

%%
%     inside=zeros(size(pt,1),1);
inside=zeros(size(pt,1),1);
num=zeros(size(pt,1),1);
for i=1:size(pt,1)  
    for j=1:size(tetra,1)  
        v1 = p(tetra(j,1),:);
        v2 = p(tetra(j,2),:);
        v3 = p(tetra(j,3),:);
        v4 = p(tetra(j,4),:);
        v=ones(4,1);
        %%
        d0=det([[v1;v2;v3;v4],v]);
        if abs(d0)>=dzero        %d0~=0                        %判断tetra的四点是否共面
            %%
            d1=det([[pt(i,:);v2;v3;v4],v]);
            d2=det([[v1;pt(i,:);v3;v4],v]);
            d3=det([[v1;v2;pt(i,:);v4],v]);
            d4=det([[v1;v2;v3;pt(i,:)],v]);
            if abs(d0-d1-d2-d3-d4)>=dzero                
                error('vic_行列式有误！');            
            end;
            d=abs([d0,d1,d2,d3,d4])<=dzero;
            dsign=sign([d0,d1,d2,d3,d4]).*(~d*1);
            %       1 1 1 1
            %面     0 1 1 1
            %棱     0 0 1 1
            %顶点   0 0 0 1
            %%
            if ~isempty(find(dsign(2:5)==-dsign(1)))        %只要有一个异号，点就在四面体外
                inside(i)=-1;
            else if ~isempty(find(dsign(2:5)==0))
                    inside(i)=0;
                else if abs(sum(dsign))==5        
                        inside(i)=1;    
                        break;
                    else error('vic_inside 无输出');
                    end
                end
            end
        else    
            error('vic_该四面体的顶点共面');    
        end
    end
    num(i)=j;
    fprintf('the %dth point\n',i)

    
    %%
%         if ~isempty(find(dsign(2:5)==-dsign(1)))        %只要有一个异号，点就在四面体外
%             inside(i)=-1;
%         else if ~isempty(find(dsign(2:5)==0))
%                 inside(i)=0;
%             else if abs(sum(dsign))==5        
%                     inside(i)=1;    
%                 else error('vic_inside 无输出');
%                 end
%             end
%         end
%%
%     if dsign(1)==1
%         inside(i,:)=dsign;
%     else
%         inside(i,:)=-dsign;            
%     end
end
    
    
    
    
    
    
    


%%
% in=p((inside==1),:);
% out=p((inside==-1),:);
% on=p((inside==0),:);
% figure;
% hold on;
% plot3(in(:,1),in(:,2),in(:,3),'.r','markersize',10);
% plot3(out(:,1),out(:,2),out(:,3),'xb','markersize',10);
% plot3(on(:,1),on(:,2),on(:,3),'ok','markersize',10);
% tetramesh(tetra,p,'edgecolor','k','FaceAlpha',0.6);
% axis([-1 5 -1 5 -1 5]);
% axis equal off;
% rotate3d on;

end

