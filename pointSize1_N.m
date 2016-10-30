function [ valid,sizeVal ] = pointSize1_N(depth,point,sizeOfDepth )
%POINTPROP Summary of this function goes here
%根据scale所生成的笛卡尔网格，矩阵大小为m1*m1*m1，计算point中各点在笛卡尔网格中的
%具体深度值
%function [ valid,size ] = pointSize1_N(depth,point,sizeOfDepth )
%input:     depth m*n*p 笛卡尔网格各点的深度值
%input:     point m2*3 空间点的坐标 
%input:     sizeOfDepth 1*3 矩阵depth的大小
%output:    valid m2*1 标记出矩阵范围外的点：如果 point的坐标在（1-n,1-m,1-p）之外，则标记为1
%output:    sizeVal m2*1
%******************************Important**********************************
%point的坐标为（x,y,z）,则对应depth中的值为depth（ceil(y),ceil(x),ceil(z)）
%******************************Important**********************************
end

