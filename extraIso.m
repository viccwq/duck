function [ ] = extraIso( tree )
%EXTRAISO Summary of this function goes here
%   Detailed explanation goes here
center=[sum(tree(:,1:2),2),sum(tree(:,3:4),2),sum(tree(:,5:6),2)]/2;
range=minmax(tree')';
xyz=center(tree(:,7)==range(2,7),:);
v=tree(tree(:,7)==range(2,7),7);

xyz=round(xyz*100);
xyz=xyz-min(min(xyz))+1;
ref=ones(size(xyz));

[xyz2,index] = sortrows(xyz,3);
ref2=ref(index);

Num=max(xyz2);
out=zeros(Num);
for i=1:Num(3)
    out(:,:,i)=full(sparse(xyz2(xyz2(:,3)==i,1),xyz2(xyz2(:,3)==i,2),ref2(xyz2(:,3)==i),Num(1),Num(2)));
end
out=(out~=0)*1+1;
[x,y,z]=ndgrid(1:Num(1),1:Num(2),1:Num(3));

p = patch(isosurface(x,y,z,out,2));
isonormals(x,y,z,out,p)
set(p,'FaceColor','red','EdgeColor','none');
daspect([1 1 1])
view(3); axis tight
camlight
lighting gouraud

end

