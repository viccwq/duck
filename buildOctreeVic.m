function [ a,b ] = buildOctreeVic( range,point )
%BUILDOCTREEVIC Summary of this function goes here
%�������㴴����άƽ��˲���
%function [ a,b ] = buildOctreeVic( range,point )
%input:    range 1*6 [xmin xmax ymin ymax zmin zmax]
%input:    point M1*4 :����Ĺ̶��� ������Ϊָ�������ֵ
%output:   a   n*8 ������Ҷ�ӵ���Ϣ[xmin��xmax��ymin��ymax��zmin��zmax��nodedepth��0/1]
%output:   b   size��point��1��*8����Ҷ�ӵ���Ϣ��[xmin��xmax��ymin��ymax��zmin��zmax��nodedepth��ID]
%debug command: 
%mex -g buildOctreeVic.cpp
%[a,b]=buildOctreeVic( [-1 1 -1 1 -1 1],node);
end

