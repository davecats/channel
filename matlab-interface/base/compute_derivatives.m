function derivatives=compute_derivatives(params)
% Compute wall-normal compact finite difference coefficients
%

y=params.y;
ny=params.ny;


M=zeros(5,5); t=zeros(5,1);

derivatives.d0=zeros(ny+3,ny+3);
derivatives.d1=derivatives.d0; 
derivatives.d2=derivatives.d0;
derivatives.d3=derivatives.d0;
derivatives.d4=derivatives.d0;


    for iY=1:ny-1
        % define helping indices
        iy=iY+2; 
        % d4
        for i=0:4; for j=0:4; M(i+1,j+1)=(y(iy-2+j)-y(iy))^(4-i); end; end
        t=t*0; t(1)=24; derivatives.d4(iy,iy-2:iy+2)=(M\t)';
        % d0
        for i=0:4; for j=0:4; M(i+1,j+1)=(5-i)*(6-i)*(7-i)*(8-i)*(y(iy-2+j)-y(iy))^(4-i); end; end
        t=t*0; for i=0:4; for j=-2:2; t(i+1)=t(i+1)+derivatives.d4(iy,iy+j)*(y(iy+j)-y(iy))^(8-i); end; end; derivatives.d0(iy,iy-2:iy+2)=(M\t)'; 
        %derivatives.d0(iy,iy-2:iy+2)=[0 0 1 0 0];
        % d3
        for i=0:4; for j=0:4; M(i+1,j+1)=(y(iy-2+j)-y(iy))^(4-i); end; end
        t=t*0; 
        for i=0:1
            for j=-2:2
                t(i+1)=t(i+1)+derivatives.d0(iy,iy+j)*(4-i)*(3-i)*(2-i)*(y(iy+j)-y(iy))^(1-i);
            end 
        end; derivatives.d3(iy,iy-2:iy+2)=(M\t)';
        % d2
        for i=0:4; for j=0:4; M(i+1,j+1)=(y(iy-2+j)-y(iy))^(4-i); end; end
        t=t*0; 
        for i=0:2 
            for j=-2:2
                t(i+1)=t(i+1)+derivatives.d0(iy,iy+j)*(4-i)*(3-i)*(y(iy+j)-y(iy))^(2-i);
            end 
        end; derivatives.d2(iy,iy-2:iy+2)=(M\t)';
        % d1
        t=t*0; 
        for i=0:3 
            for j=-2:2
                t(i+1)=t(i+1)+derivatives.d0(iy,iy+j)*(4-i)*(y(iy+j)-y(iy))^(3-i);
            end 
        end; derivatives.d1(iy,iy-2:iy+2)=(M\t)';
    end
    % bottom wall
    iY=0; iy=iY+2;
    for i=0:4; for j=0:4; M(i+1,j+1)=(y(iy-2+j+1)-y(iy))^(4-i); end; end
    t=t*0; t(4)=1; derivatives.d1(iy,1:5)=(M\t)';
    t=t*0; t(3)=2; derivatives.d2(iy,1:5)=(M\t)';
    t=t*0; t(2)=6; derivatives.d3(iy,1:5)=(M\t)';
    derivatives.d0(iy,iy)=1;
    % bottom ghost
    iY=-1; iy=iY+2;
    for i=0:4; for j=0:4; M(i+1,j+1)=(y(iy-2+j+2)-y(iy))^(4-i); end; end
    t=t*0; t(4)=1; derivatives.d1(iy,1:5)=(M\t)';
    t=t*0; t(3)=2; derivatives.d2(iy,1:5)=(M\t)';
    t=t*0; t(2)=6; derivatives.d3(iy,1:5)=(M\t)';
    derivatives.d0(iy,iy)=1;
    % top wall
    iY=ny; iy=iY+2;
    for i=0:4; for j=0:4; M(i+1,j+1)=(y(iy-2+j-1)-y(iy))^(4-i); end; end
    t=t*0; t(4)=1; derivatives.d1(iy,end-4:end)=(M\t)';
    t=t*0; t(3)=2; derivatives.d2(iy,end-4:end)=(M\t)';
    t=t*0; t(2)=6; derivatives.d3(iy,end-4:end)=(M\t)';
    derivatives.d0(iy,iy)=1;
    % top ghost
    iY=ny+1; iy=iY+2;
    for i=0:4; for j=0:4; M(i+1,j+1)=(y(iy-2+j-2)-y(iy))^(4-i); end; end
    t=t*0; t(4)=1; derivatives.d1(iy,end-4:end)=(M\t)';
    t=t*0; t(3)=2; derivatives.d2(iy,end-4:end)=(M\t)';
    t=t*0; t(2)=6; derivatives.d3(iy,end-4:end)=(M\t)';
    derivatives.d0(iy,iy)=1;
