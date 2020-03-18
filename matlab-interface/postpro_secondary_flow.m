%% Cleanup environment
%  ------------------------------
clear variables
close all
addpath('./base/');

%% Plot setting
%  ------------------------------
plots.lw=2;

%% Init simulation parameters
%  ------------------------------
dns = init_dns('../dns.in',1);
der = compute_derivatives(dns);

%% Open file with averages
%  ------------------------------
data=open_binary('../cpl-interface/post.bin');

%% Body contour
%  ------------------------------
f=fopen('../ibm.bin','r');
field = permute(reshape(fread(f,'double'),dns.sized(3:-1:1)),[3,2,1]);
fclose(f);
B = bwboundaries(field(:,:,1)<1,4);

%% Phase average
%  ------------------------------
dns.Lz = 2*pi/dns.beta0;
ibm.h = 0.02352;
ibm.s = 0.04912;
ibm.base = ibm.s*148.188*2/614;
ibm.Nrib = round(dns.Lz/ibm.s);
ibm.NN = 2*ceil(dns.nzd/ibm.Nrib);
ibm.zz = (0:ibm.NN-1)*ibm.s/ibm.NN;
for j=1:numel(data)
    sz=size(data{j}.data);
    DATAi=fft(data{j}.data,[],numel(sz));
    kz=repmat(reshape([(0:floor(dns.nzd/2)-1) 0 (-floor(dns.nzd/2)+1:-1)]*dns.beta0,[ones(1,numel(sz)-1) sz(end)]),[sz(1:end-1) 1]);
    DATAo=DATAi;
    for iRib = 1:ibm.Nrib-1
        DATAo=DATAo+DATAi.*exp(-1j*ibm.s*iRib*kz);
    end
    data{j}.avgdata=zeros([sz(1:end-1) ibm.NN]);
    for iz=1:ibm.NN
        idx=repmat({':'}, 1, numel(sz)-1);
        data{j}.avgdata(idx{:},iz)=sum(DATAo(idx{:},:).*exp(kz(idx{:},:)*1j*ibm.zz(iz)),numel(sz))/dns.nzd^(1.5);
    end
end

%% One dimensional analysis
%  ------------------------------
%  ------------------------------

%% Actual channel height
%  ------------------------------
dns.H = dns.y(end-1)-dns.y(2);
dns.Hmelt = dns.H*(1-ibm.base/ibm.s) + (dns.H-ibm.h/2)*ibm.base/ibm.s;
dns.yH = (dns.y - dns.y(2))/dns.H;

%% Ubulk
%  ------------------------------
stats.Ub = trapz(dns.y,(mean(data{1}.data(1,:,:),3)))/dns.Hmelt;
stats.Reb = stats.Ub*dns.H/dns.nu;

%% Total stresses
%  ------------------------------

% Reynolds shear stress
stats.profiles.uv = (mean(data{2}.data(4,:,:),3));
% d<U>/dy
stats.profiles.Uy = (mean(data{3}.data(2,:,:),3)); % XXX do we need a division by H here?
% Immerse boundary force (needed for intrinsic averaging)
%D = der.d1; 
%D(end-1,:)=0; D(end-1,end-1)=1;
%stats.profiles.intfx = (D\(der.d0*mean(data{4}.data(1,:,:),3)'))';
stats.profiles.intfx = fliplr(cumtrapz(dns.yH(end:-1:1),mean(data{4}.data(1,end:-1:1,:),3)));
stats.profiles.totshear = -stats.profiles.uv+dns.nu*stats.profiles.Uy+stats.profiles.intfx;

figure(99)
hold on; box on
plot(dns.y,dns.nu*stats.profiles.Uy,'-','Linewidth', plots.lw)
plot(dns.y,-stats.profiles.uv, '-','Linewidth', plots.lw)
plot(dns.y,stats.profiles.intfx, '-','Linewidth', plots.lw)
plot(dns.y,stats.profiles.totshear, '-','Linewidth', plots.lw, 'Color', 'k')
plot(dns.y, (1-dns.y))
plot([0 0], [0 1], 'k--')
xlim([dns.y(1) 1])
ylim([0 1])
legend('d<U>/dy', '<uv>', '\int f_x', 'tot' , 'linear')

%% Retau (I am assuming here that the equivalent wall shear stress, so utau is unity)
%  ------------------------------
stats.retau = dns.Hmelt/dns.nu;

return

%% Two dimensional analysis
%  ------------------------------
%  ------------------------------

%% Viscous stress along surface (with interpolation where the body should be
%  ------------------------------

% Surface and normal definition
contour.N=601;
contour.z = linspace(0,ibm.s,contour.N);
contour.y = dns.y(1)+(contour.z<0.5*(ibm.s+ibm.base)).*(contour.z>0.5*(ibm.s-ibm.base)).*(ibm.h-ibm.h*abs(contour.z-ibm.s/2)/ibm.base*2);
contour.n = zeros(3,numel(contour.z));
contour.n(2,:) = 1-(contour.z<0.5*(ibm.s+ibm.base)).*(contour.z>0.5*(ibm.s-ibm.base)).*(1-0.5*ibm.base/(sqrt(ibm.h^2+(ibm.base/2)^2)));
contour.n(2,floor(contour.N/2)+1)=1;
contour.n(3,:) = (contour.z<0.5*(ibm.s+ibm.base)).*(contour.z>0.5*(ibm.s-ibm.base)).*(ibm.h/(sqrt(ibm.h^2+(ibm.base/2)^2))).*sign(contour.z-ibm.s/2);
contour.s = cumsum(sqrt([0 diff(contour.z)].^2+[0 diff(contour.y)].^2));

% Interpolate data
[ZZ,YY] = ndgrid(ibm.zz,dns.y);
FsigmaUy = griddedInterpolant(ZZ, YY, real(squeeze(data{3}.avgdata(2,:,:)))'); 
FsigmaUz = griddedInterpolant(ZZ, YY, real(squeeze(data{3}.avgdata(3,:,:)))'); 
FsigmaWy = griddedInterpolant(ZZ, YY, real(squeeze(data{3}.avgdata(8,:,:)))'); 
FsigmaWz = griddedInterpolant(ZZ, YY, real(squeeze(data{3}.avgdata(9,:,:)))'); 
stats.slices.sigmaUy=FsigmaUy(contour.z,contour.y);
stats.slices.sigmaUz=FsigmaUz(contour.z,contour.y);
stats.slices.sigmaWy=FsigmaWy(contour.z,contour.y);
stats.slices.sigmaWz=FsigmaWz(contour.z,contour.y);
stats.slices.sigma_x=stats.slices.sigmaUy.*contour.n(2,:)+stats.slices.sigmaUz.*contour.n(3,:);

tau_w = dns.nu*trapz(contour.s/ibm.s,stats.slices.sigma_x)

figure(199)
hold on
plot3(contour.z,contour.y,(stats.slices.sigma_x)*dns.nu)
plot3(contour.z,contour.y,(stats.slices.sigmaUy)*dns.nu.*contour.n(2,:))
plot3(contour.z,contour.y,(stats.slices.sigmaUz)*dns.nu.*contour.n(3,:))

%figure(); hold on; plot(ibm.zz,(((squeeze(real(data{3}.avgdata(3,80,:))))))); plot(ibm.zz,smooth(((squeeze(real(data{3}.avgdata(3,80,:))))),6))


%% Viscous stress along surface (wall located at maximum shear stress position)
%  ------------------------------

mid = stats.slices.sigmaUy(floor(end/2)+1);

% Surface and normal definition
contour.y=unique(dns.y(dns.y<=(ibm.h+dns.y(2)))); contour.y=contour.y(3:end);
stats.slices.sigmaUy=zeros(size(contour.y));
Uyl=zeros(size(contour.y)); Uyr=zeros(size(contour.y));
Uzl=zeros(size(contour.y)); Uzr=zeros(size(contour.y));
for i=1:numel(contour.y)
    window_size=1;
    % Uzl
    smoothed = smooth(((squeeze(real(data{3}.avgdata(3,i+2,:))))),window_size);
    [m,im] = sort(abs(smoothed),'descend');
    sm = sign(smoothed(im));
    Uzl(i) = sm(1)*m(1); iim=find(sm>0,1); Uzr(i) = sm(iim)*m(iim);
    % Uzl
    smoothed = smooth(((squeeze(real(data{3}.avgdata(2,i+2,:))))),window_size);
    Uyl(i) = smoothed(im(1)); Uyr(i) = smoothed(im(iim));
    % 
    zl(i) = ibm.zz(im(1)); zr(i) = ibm.zz(im(iim));
end
i0 = find(ibm.zz/ibm.s*2==1);
[~, iy0] = max(abs(real(data{3}.avgdata(2,:,i0))));

il = find(ibm.zz<=(ibm.s/2 - ibm.base/2));
ir = find(ibm.zz>=(ibm.s/2 + ibm.base/2));
contour.z = [ ibm.zz(il) zl ibm.zz(i0) fliplr(zr) ibm.zz(ir)];
contour.y = [ dns.y(2)*ones(1,numel(il)) contour.y ibm.h fliplr(contour.y) dns.y(2)*ones(1,numel(ir)) ];
stats.slices.sigmaUy = [ squeeze(real(data{3}.avgdata(2,2,il)))' Uyl real(data{3}.avgdata(2,iy0,i0))  fliplr(Uyr) squeeze(real(data{3}.avgdata(2,2,ir)))' ];
stats.slices.sigmaUz = [ squeeze(real(data{3}.avgdata(3,2,il)))' Uzl real(data{3}.avgdata(3,iy0,i0)) fliplr(Uzr) squeeze(real(data{3}.avgdata(3,2,ir)))' ];
ny0= 0.5*ibm.base/(sqrt(ibm.h^2+(ibm.base/2)^2)); contour.ny = [ones(1,numel(il))  ny0*[ones(1,numel(zl)) 1/ny0 ones(1,numel(zl))] ones(1,numel(ir))];
nz0= (ibm.h/(sqrt(ibm.h^2+(ibm.base/2)^2)));      contour.nz = [zeros(1,numel(il)) nz0*[ones(1,numel(zl)) 0 ones(1,numel(zl))] zeros(1,numel(ir))];
contour.nz = contour.nz.*sign(contour.z-ibm.s/2);
contour.s = cumsum(sqrt([0 diff(contour.z)].^2+[0 diff(contour.y)].^2));
stats.slices.sigma_x = (stats.slices.sigmaUy.*contour.ny+stats.slices.sigmaUz.*contour.nz);

tau_w = dns.nu*trapz(contour.s/ibm.s,stats.slices.sigma_x)

figure(199)
hold on
plot3(contour.z,contour.y,(stats.slices.sigma_x)*dns.nu,'.')
plot3(contour.z,contour.y,(stats.slices.sigmaUy)*dns.nu.*contour.ny)
plot3(contour.z,contour.y,(stats.slices.sigmaUz)*dns.nu.*contour.nz)

%figure(); hold on; i=3; iy=126; plot(ibm.zz,(((squeeze(real(data{3}.avgdata(i,iy,:))))))); plot(ibm.zz,smooth(((squeeze(real(data{3}.avgdata(i,iy,:))))),8))


%% Plot (z-y) slice 
figure(1)
hold on 
box on
%surf(dns.z,dns.y, squeeze(sum(data{1}.data(1:3,:,:).^2,1))); shading interp; 
surf(dns.z,dns.y, ((squeeze(real(data{1}.data(1,:,:)))))); shading interp; 
view([0 90])
colorbar
set(gca(),'Layer','top')
%contour(dns.z,dns.y,field(:,:,1),[1 1],'Linewidth',2,'Color','w')
plot3(dns.z(B{1}(:,2)),dns.y(B{1}(:,1)),1000*ones(size(B{1}(:,1))), 'w','linewidth',2)
% streamslice(dns.z,dns.y, squeeze(data{1}.data(3,:,:)),  squeeze(data{1}.data(2,:,:)),80);

%% Plot (z-y) conditional average slice
figure(10)
set(gcf(),'Unit','Centimeters','Position',[0 0 6 25])
hold on 
box on
%surf(dns.z,dns.y, squeeze(sum(data{1}.data(1:3,:,:).^2,1))); shading interp; 
surf(ibm.zz,dns.y, ((squeeze(real(data{1}.avgdata(2,:,:)))))); shading interp; 
ZZ = dns.z(B{1}(:,2)); YY = dns.y(B{1}(:,1));
Z = ZZ(ZZ<ibm.s & ZZ>0 & YY<1);
Y = YY(ZZ<ibm.s & ZZ>0 & YY<1);
plot3(Z,Y,1000*ones(size(Z)), 'w--','linewidth',2)
ylim([dns.y(2) 1])
xlim([0 ibm.s])
view([0 90])
colorbar
set(gca(),'Layer','top')
%contour(dns.z,dns.y,field(:,:,1),[1 1],'Linewidth',2,'Color','w')
%plot3(dns.z(B{1}(:,2)),dns.y(B{1}(:,1)),1000*ones(size(B{1}(:,1))), 'w','linewidth',2)
% streamslice(dns.z,dns.y, squeeze(data{1}.data(3,:,:)),  squeeze(data{1}.data(2,:,:)),80);

%%
figure(2)
hold on
plot((dns.y)*1000,(mean(data{1}.data(1,:,:),3)))
plot((dns.y)*1000,(1/0.39)*log((dns.y)*1000)+4.5)
set(gca(),'Xscale','log')
% figure(2)
% hold on
% % dwdz
% SUB = [9*ones(numel(B{1}(:,1)),1) B{1}(:,1)  B{1}(:,2)];
% IND = sub2ind(size(data{3}.data), SUB(:,1), SUB(:,2), SUB(:,3));
% plot3(dns.z(B{1}(:,2)),dns.y(B{1}(:,1)),data{3}.data(IND))
% % dudy
% SUB = [2*ones(numel(B{1}(:,1)),1) B{1}(:,1)  B{1}(:,2)];
% IND = sub2ind(size(data{3}.data), SUB(:,1), SUB(:,2), SUB(:,3));
% plot3(dns.z(B{1}(:,2)),dns.y(B{1}(:,1)),data{3}.data(IND))
% 
% figure(3)
% dudy=zeros(dns.ny/2+1,dns.nzd);
% deltay=zeros(1,dns.nzd);
% for iz=1:dns.nzd
%     iy=find(field(:,iz,1)<=0.5,1);
%     dudy(:,iz)=abs(data{3}.data(2,iy,iz))+abs(data{3}.data(9,iy,iz));
%     deltay(iz)=dns.y(iy);
% end
% utau=sqrt(abs(dudy*dns.nu));
% hold on
% set(gca(),'Xscale','log')
% for iz=1:dns.nzd
%   if utau(1,iz)>1
%     semilogx((dns.y(2:floor(dns.ny/2)+2))/(1-deltay(iz))-deltay(iz),(squeeze(data{1}.data(1,2:floor(dns.ny/2)+2,iz))-data{1}.data(1,floor(dns.ny/2)+2,iz))./utau(:,iz)')
%   end
% end
% Y=dns.y(2:floor(dns.ny/2)+2);
% plot(Y, (1/0.395)*log(Y))
% 
% 
% 
% figure(4)
% dudy=zeros(dns.ny/2+1,dns.nzd);
% deltay=zeros(1,dns.nzd);
% for iz=1:dns.nzd
%     iy=find(field(:,iz,1)<=0.5,1);
%     dudy(:,iz)=abs(data{3}.data(2,iy,iz))+abs(data{3}.data(9,iy,iz));
%     deltay(iz)=dns.y(iy);
% end
% utau=sqrt(abs(dudy*dns.nu));
% hold on
% set(gca(),'Xscale','log')
% for iz=1:dns.nzd
%   if utau(1,iz)>1
%     if deltay(iz)>0
%         color='r';
%     else
%         color='b';
%     end
%     semilogx((dns.y(2:floor(dns.ny/2)+2)-deltay(iz))*utau(1,iz)/dns.nu,(squeeze(data{1}.data(1,2:floor(dns.ny/2)+2,iz)))./utau(:,iz)',color)
%   end
% end
% xlim([1 300])
% Y=dns.y(2:floor(dns.ny/2)+2);
% plot(Y*200, (1/0.395)*log(Y*200),'k')
