function [ valid,sizeVal ] = pointSize1_N(depth,point,sizeOfDepth )
%POINTPROP Summary of this function goes here
%����scale�����ɵĵѿ������񣬾����СΪm1*m1*m1������point�и����ڵѿ��������е�
%�������ֵ
%function [ valid,size ] = pointSize1_N(depth,point,sizeOfDepth )
%input:     depth m*n*p �ѿ��������������ֵ
%input:     point m2*3 �ռ������� 
%input:     sizeOfDepth 1*3 ����depth�Ĵ�С
%output:    valid m2*1 ��ǳ�����Χ��ĵ㣺��� point�������ڣ�1-n,1-m,1-p��֮�⣬����Ϊ1
%output:    sizeVal m2*1
%******************************Important**********************************
%point������Ϊ��x,y,z��,���Ӧdepth�е�ֵΪdepth��ceil(y),ceil(x),ceil(z)��
%******************************Important**********************************
end

