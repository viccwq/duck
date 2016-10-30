function showSta(p,tetra)
%SHOWSTA Summary of this function goes here
%������ϵ��Ϊ�����弸�����ġ����������֮��ľ����������뾶�ı�ֵ����
%ֵ�� 0 �� 1 ֮��
%��ʾ��Ԫ������
out=zeros(size(tetra,1),1);
count=zeros(1,100);
for i=1:size(tetra,1)
    point=[p(tetra(i,1),:);p(tetra(i,2),:);p(tetra(i,3),:);p(tetra(i,4),:)];
    center=sum(point)/4;
    val=fcircumsphere(point);
    out(i)=sum((val(1:3)-center).^2,2)/val(4);
    index=round(out(i)*100);
    if(index==0)
        index=1;
    end
    count(index)=count(index)+1;
end
figure(2);
bar(1:100,count,1,'c');
title('���ġ����������֮��ľ����������뾶�ı�ֵ');
% disp('press any key to continue');
% pause;
end

