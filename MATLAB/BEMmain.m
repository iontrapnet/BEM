clear;%改变points前一定要删除或转移目录下Ve.mat
Vf=[1 0 0 0 0 1];%每个电极电势
M=csvread('4rod.csv');
noe=M(1,2);
p=size(M,1);
%载入计算点
xyzmin=M(p-2,:);
xyzmax=M(p-1,:);
npoint=M(p,:);
if npoint(1)==1
    xp=xyzmin(1);
else
    xp=xyzmin(1):((xyzmax(1)-xyzmin(1))/npoint(1)):xyzmax(1);
end
if npoint(1)==1
    yp=xyzmin(2);
else
    yp=xyzmin(2):((xyzmax(2)-xyzmin(2))/npoint(2)):xyzmax(2);
end
if npoint(3)==1
    zp=xyzmin(3);
else
    zp=xyzmin(3):((xyzmax(3)-xyzmin(3))/npoint(3)):xyzmax(3);
end
points=zeros(size(xp,2)*size(yp,2)*size(zp,2),3);
ijk=1;
for i=1:size(xp,2)
    for j=1:size(yp,2)
        for k=1:size(zp,2)
            points(ijk,:,:)=[xp(i),yp(j),zp(k)];
            ijk=ijk+1;
        end
    end
end
%查找是否计算过
if exist('Ve.mat','file')==2
    load('Ve.mat');
else
    BEMsub(points);
end
V=zeros(size(points,1),1);
for i=1:noe
    V=Ve(:,i).*Vf(i)+V;
end