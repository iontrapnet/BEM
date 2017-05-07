function BEMsub(points)
%改变原始数据点前一定要删除或转移目录下charge.mat
M=csvread('4rod.csv');
non=M(1,1);
noe=M(1,2);
nodes=M(2:non+1,:);
not=zeros(noe,1);
p=size(M,1);
V=zeros(p-non-noe-4,noe);
Element=zeros(p-non-noe-4,3);
% xyzmin=M(p-2,:);
% xyzmax=M(p-1,:);
% npoint=M(p,:);
%Vel=[1 1 1 1 1 1];%每个电极电势
% if npoint(1)==1
%     xp=xyzmin(1);
% else
%     xp=xyzmin(1):((xyzmax(1)-xyzmin(1))/npoint(1)):xyzmax(1);
% end
% if npoint(1)==1
%     yp=xyzmin(2);
% else
%     yp=xyzmin(2):((xyzmax(2)-xyzmin(2))/npoint(2)):xyzmax(2);
% end
% if npoint(3)==1
%     zp=xyzmin(3);
% else
%     zp=xyzmin(3):((xyzmax(3)-xyzmin(3))/npoint(3)):xyzmax(3);
% end
%points=zeros(size(xp,2)*size(yp,2)*size(zp,2),3);
% ijk=1;
% for i=1:size(xp,2)
%     for j=1:size(yp,2)
%         for k=1:size(zp,2)
%             points(ijk,:,:)=[xp(i),yp(j),zp(k)];
%             ijk=ijk+1;
%         end
%     end
% end

%按电极给出原始点和不同的初始边界
a=non+3;
av=1;
bv=0;
for i=1:noe
    not(i)=M(non+2+sum(not),1)+1;
    b=non+1+sum(not);
    mn=size(M(a:b,:),1);
    bv=mn+bv;
    V(av:bv,i)=ones(mn,1);
    Element(av:bv,:)=M(a:b,:);
    av=bv+1;
    a=b+2;
end

%计算charge
len=size(Element,1);
x=nodes(:,1);
y=nodes(:,2);
z=nodes(:,3);
node1=Element(:,1);
node2=Element(:,2);
node3=Element(:,3);
xn=zeros(1,len);
yn=zeros(1,len);
zn=zeros(1,len);
triangles=zeros(3,3,len);
for j=1:len
    xn(j)=(x(node1(j))+x(node2(j))+x(node3(j)))/3;
    yn(j)=(y(node1(j))+y(node2(j))+y(node3(j)))/3;
    zn(j)=(z(node1(j))+z(node2(j))+z(node3(j)))/3;
    triangles(:,:,j)=[[x(node1(j)),y(node1(j)),z(node1(j))];[x(node2(j)),y(node2(j)),z(node2(j))];[x(node3(j)),y(node3(j)),z(node3(j))]];
end
%查找是否已经计算过
if exist('charge.mat','file')==2
    load('charge.mat');
else
    A=[xn;yn;zn]';
    alpha=Malpha(triangles,A,len);
    charge=zeros(len,noe);
    for i=1:noe
        chargep=alpha\V(:,i);%不含4pi
        charge(:,i)=chargep.*4.*pi;
    end
    save('charge.mat','charge');
end
Ve=zeros(size(points,1),noe);
for i=1:noe
    Ve(:,i) = potential(triangles, charge(:,i), points);
end
save('Ve.mat','Ve');
end