function [cb,triangles] = ChargeBasis(path)
    file=[path '.mat'];
    if exist(file,'file')==2
        load(file);
        return;
    end
    
    M=csvread([path '.csv']);
    non=M(1,1);
    noe=M(1,2);
    nodes=M(2:non+1,:);
    not=zeros(noe,1);
    p=size(M,1);
    V=zeros(p-non-noe-4,noe);
    Element=zeros(p-non-noe-4,3);
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

    A=[xn;yn;zn]';
    alpha=Malpha(triangles,A,len);
    cb=zeros(len,noe);
    for i=1:noe
        chargep=alpha\V(:,i);
        cb(:,i)=chargep.*4.*pi;
    end
    
    save(file,'cb','triangles');
end