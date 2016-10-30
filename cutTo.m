function vec_new=cutTo(ori_node,end_node,vec,depth,Dep)
%CUTTO Summary of this function goes here
%   Detailed explanation goes here
if Dep>=0
    error('Dep����С��0');
end
count=9;
vec_new=zeros(size(vec));

temp=repmat([1:count]',1,3)/(count+1);
for i=1:size(ori_node,1)
    vec_temp=repmat(vec(i,:),count,1).*temp;
    point=repmat(ori_node(i,:),count,1)+vec_temp;
    [ outside,sizeVal ] = pointSize1_N(depth,point,size(depth));
    j=sum(sizeVal>=(Dep+max(abs(diff(sizeVal)))));
    vec_new(i,:)=vec_temp(j,:);
end

end

