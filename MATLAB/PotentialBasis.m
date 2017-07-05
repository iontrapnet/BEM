function pb=PotentialBasis(path, xr, yr, zr)
    global triangles
    file=[path '-' DataHash([xr yr zr]) '.mat'];
    if exist(file,'file')==2
        load(file);
        return;
    end
    cb = ChargeBasis(path);
    [y,z,x]=meshgrid(linspace(xr(1),xr(2),xr(3)+1),linspace(yr(1),yr(2),yr(3)+1),linspace(zr(1),zr(2),zr(3)+1));
    nop=[numel(x) 1];
    points=[reshape(x,nop) reshape(y,nop) reshape(z,nop)];
    noe=size(cb,2);
    
    % single machine
    pb=zeros(4,nop(1),noe);
    
    % parallel
    %part=max([xr(3) yr(3) zr(3)])+1;
    %points=reshape(points,[part,nop(1)/part,3]);
    %pb=zeros(4,part,nop(1)/part,noe);
    %parfor i=1:part
    
        ve=eye(noe);
        %parfor k=1:noe
        for k=1:noe
            
            % single machine
            [pot,fx,fy,fz]=Potential(triangles,cb*ve(:,k),points);
            pb(:,:,k)=[pot;fx;fy;fz];
            
            % parallel
            %[pot,fx,fy,fz]=Potential(triangles,cb*ve(:,k),shiftdim(points(i,:,:)));
            %pb(:,i,:,k)=[pot;fx;fy;fz];
        end
    
    % parallel
    %end
    %pb=reshape(pb,[4,nop(1),noe]);
    
    save(file,'pb');
end