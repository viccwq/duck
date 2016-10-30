function [ output ] = movePoint(depth,ori_node,end_node,norm,vec)
%MOVEPOINT Summary of this function goes here
%   Detailed explanation goes here
if dot(norm,vec)>=0
    if dot(norm,vec)>0
        norm=-1*norm;
    else
        pause;
    end    
end


a=repmat(end_node,8,1)+repmat([1;2;3;4;5;6;7;8]*1,1,3).*repmat(norm,8,1);
[ valid,sizeVal] = pointSize1_N(depth,a,size(depth));
test=sizeVal(1:end-1)>sizeVal(2:end);
if sum(test)>=7
     norm=-1*norm;
end
k=1;
p=end_node+k*norm;
[ valid,sizeVal ] = pointSize1_N(depth,p,size(depth));
count=0;
while abs(sizeVal)>0.1
    if sizeVal>0
        k=k/2;
        p=p-k*norm;
    else
        if k~=1
            k=k/2;
        end
        p=p+k*norm;       
    end
    [ valid,sizeVal ] = pointSize1_N(depth,p,size(depth));
%     if count==40
%         [MCon,NCon,PCon]=size(depth);
%         %x方向
%         if ori_point(1)<=1
%             ori_point(1)=3;
%         end
%         if ori_point(1)>=NCon                    
%             ori_point(1)=NCon-2;
%         end
%          %y方向
%        if ori_point(2)<=1;
%             ori_point(2)=3;
%         end
%         if ori_point(2)>=MCon                    
%             ori_point(2)=MCon-2;
%         end
%         %z方向
%         if ori_point(3)<=1;
%             ori_point(3)=3;
%         end
%         if ori_point(3)>=PCon                    
%             ori_point(3)=PCon-2;
%         end
%         p=ori_point;
%         break;
%     end
%     count=count+1;
end
output=p;
end

