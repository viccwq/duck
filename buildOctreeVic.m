function [ a,b ] = buildOctreeVic( range,point )
%BUILDOCTREEVIC Summary of this function goes here
%由所给点创建三维平衡八叉树
%function [ a,b ] = buildOctreeVic( range,point )
%input:    range 1*6 [xmin xmax ymin ymax zmin zmax]
%input:    point M1*4 :输入的固定点 第四列为指定的深度值
%output:   a   n*8 点所在叶子的信息[xmin，xmax，ymin，ymax，zmin，zmax，nodedepth，0/1]
%output:   b   size（point，1）*8所有叶子的信息：[xmin，xmax，ymin，ymax，zmin，zmax，nodedepth，ID]
%debug command: 
%mex -g buildOctreeVic.cpp
%[a,b]=buildOctreeVic( [-1 1 -1 1 -1 1],node);
end

