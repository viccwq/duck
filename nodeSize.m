
function [eledist]=nodeSize(node,ele)
%%
%[eledist]=nodeSize(node,ele)
%获取每个节点的实际尺寸
%输入：node        n*3
%输入：ele         m*3
%输出：eledist     n*1

%%
%1-2
d1=sqrt((node(ele(:,1),1)-node(ele(:,2),1)).^2+...
        (node(ele(:,1),2)-node(ele(:,2),2)).^2+...
        (node(ele(:,1),3)-node(ele(:,2),3)).^2);
%1-3
d2=sqrt((node(ele(:,1),1)-node(ele(:,3),1)).^2+...
        (node(ele(:,1),2)-node(ele(:,3),2)).^2+...
        (node(ele(:,1),3)-node(ele(:,3),3)).^2);
%2-3
d3=sqrt((node(ele(:,2),1)-node(ele(:,3),1)).^2+...
        (node(ele(:,2),2)-node(ele(:,3),2)).^2+...
        (node(ele(:,2),3)-node(ele(:,3),3)).^2);    

eleindex=ele(:,[1,2,1,3,2,3]);
k=sparse(eleindex,ones(size(eleindex)),ones(size(eleindex)));
k=full(k);
eledist=sparse(eleindex,ones(size(eleindex)),[d1,d1,d2,d2,d3,d3]);
eledist=full(eledist);
eledist=eledist./k;

return
end
