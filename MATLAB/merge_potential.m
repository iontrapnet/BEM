cd D:\GitHub\iontrapnet\BEMParallel\MATLAB;

path='..\\Model\\4rod\\167634622912717531';
xr=[-0.005 0.005 100];
yr=[-0.005 0.005 100];
zr=[2.095 2.105 100];
%FieldInit(path,xr,yr,zr);

div=25;
nx=xr(3)+1;
ny=yr(3)+1;
nz=zr(3)+1;
noe=6;
pbs=zeros(4,nx*ny*nz,noe);
gap=(zr(2)-zr(1))/zr(3);
for part=1:(zr(3)/div)
    file = [path '-' DataHash([xr yr [zr(1)+div*(part-1)*gap zr(1)+div*part*gap div]]) '.mat'];
    if exist(file,'file')==2
        %whos('-file',file)
        load(file);
        for i=0:nx-1
            for j=0:ny-1
                for k=0:div 
                    pbs(:,i*ny*nz+j*nz+div*(part-1)+k+1,:)=pb(:,i*ny*(div+1)+j*(div+1)+k+1,:);
                end
            end
        end
    end
end

file = [path '-' DataHash([xr yr zr]) '.mat'];
load(file);