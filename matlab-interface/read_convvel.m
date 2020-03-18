%% Cleanup environment
%  ------------------------------
clear variables
close all
addpath('./base/');

%% Init simulation parameters
%  ------------------------------
dns = init_dns('../dns.in');
nfmin=85;
nfmax=120;

%% Open file 
%
sz = [dns.nzd,dns.nx+1,dns.ny+3,3];
data=zeros(sz);
for i=nfmin:nfmax
    fname=strcat('../Convvel.cart.',num2str(i),'.out')
    f = fopen(fname,'r');
    data = data + reshape(fread(f,'double'), sz);
    fclose(f);
end
data=data/(nfmax-nfmin+1);

%% Make further averages
%  ------------------------------
Nphases = 4; 
data(:,:,:)=0.5*(data(:,:,:)+data(:,:,end:-1:1)); 
%
for i=1:Nphases-1
    data=0.5*(data+circshift(data,-i*dns.nzd/Nphases,1));
end
data=circshift(data,-dns.nzd/Nphases/8,1);
data((i-1)*dns.nzd/Nphases+1:i*dns.nzd/Nphases,:,:)=0.5*(data((i-1)*dns.nzd/Nphases+1:i*dns.nzd/Nphases,:,:)+data(i*dns.nzd/Nphases:-1:(i-1)*dns.nzd/Nphases+1,:,:));
data=circshift(data,dns.nzd/Nphases/8,1);


%% Body contour
field = zeros(dns.sized(1:2),'double');
%
ibm.h = dns.y(25);
ibm.N = 4;
for ix=1:1 %dns.sized(3)
    %disp([num2str(ix) ' of ' num2str(dns.sized(3))]);
    v = zeros(dns.sized(1:2)); N=ibm.N;
    for iz=1:dns.sized(2)
        for iy=1:dns.sized(1)
            v(iy,iz) = (mod(iz+dns.nzd/N/2-1,dns.nzd/N)<48)*( (dns.y(iy)<=ibm.h)*(dns.y(iy)<1) + ((dns.ymax-dns.y(iy))<ibm.h)*(dns.y(iy)>=1) );
        end
    end
    field(:,:,ix)=v;
end
B = bwboundaries(field(:,:)<0.5,4);

%% Plot
ix=50;
iV=2;
figure(ix)
set(gcf(),'Units','centimeters','Position',[0 0 1 1.5]*7);
hold on
box on
pcolor(dns.z,dns.y,-((squeeze(data(:,ix,:,iV))))'); shading interp
plot(dns.z(B{1}(:,2)),dns.y(B{1}(:,1)),'k', 'linewidth',2)
set(gca(),'Clim',[8 25])
xlim([0.5+0.25/2-0.5 0.5+0.25/2+0.5])
ylim([0 2])
set(gca(),'Layer','top','Linewidth',2)
colorbar();