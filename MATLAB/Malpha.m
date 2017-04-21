function alpha=Malpha(triangles,A,len)
alpha=zeros(len);
for j=1:len
    alpha(:,j)=int_green3d_tri(A,triangles(:,:,j));
end%alpha矩阵没有乘1/4/pi
%normald=alpha\V;%不含4pi
%normaldp=normald.*4.*pi;%标准表面电荷
end