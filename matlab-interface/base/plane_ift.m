function dataOUT=plane_ift(dataIN,dataOUT,xdim,zdim)
% Backward fourier transform

sin = size(dataIN);  
ndims = length(sin); idx=repmat({':'}, 1, ndims); dims=[xdim zdim];
nz=sin(dims(2));   nx=sin(dims(1));   nz=(nz-1)/2; nx=nx-1;
sout= size(dataOUT); nzd=sout(dims(2)); nxd=sout(dims(1)); nxd=nxd/2;

% Swap the two half of the vector from (-n..n) --> (0..n, -n..-1)
idxOUT = idx; idxIN=idx; 
idxOUT{dims(2)}=[nzd+(-nz:-1)+1 1+(0:nz)];  idxIN{dims(2)}=[(-nz:-1)+nz+1 (0:nz)+nz+1];
idxOUT{dims(1)}=(0:nx)+1;                   idxIN{dims(1)}=(0:nx)+1;
dataOUT(idxOUT{:}) = dataIN(idxIN{:});

% Backward fourtier transform
idxOUT = idx; idxIN=idx; 
idxOUT{dims(1)}=(0:nx)+1;        idxIN{dims(1)}=(0:nx)+1;
dataOUT(idxOUT{:}) = ifft(dataOUT(idxIN{:}),[],dims(2));

% 1D transform in the axial direction
dataOUT = ifft(dataOUT,[],dims(1),'symmetric')*(2*nzd*nxd);


end