close all

%% INPUT
filename = './gke.bin';
ny=500;       nx=20;          nz=82;
alfa0=15.708; beta0=6.28319;  a=1.66;
Re=1000;
premultiply=true;

%% Compute coordinates
y=tanh(a*(2*(-1:ny+1)/ny-1))/tanh(a)+1; y=y(1:ny/2+2);
kz=(0:nz)*beta0; lambdaz=2*pi./kz;
f=fopen('y.bin'); y=fread(f,'double');

%% Load file
f=fopen('uiuj_spectra.bin','r');
mean=fread(f,[7,ny+3],'double');
uiuj.spectra=reshape(fread(f,[(2*nz+1)*10*6,ny+3],'double'),[6,10,2*nz+1,ny+3]);
uiuj.profile=reshape(fread(f,[10*6,ny+3],'double'),[6,10,ny+3]);
uiuj.convs=reshape(fread(f,[4*(2*nz+1)*(nz+1),ny+3],'double'),[4,2*nz+1,nz+1,ny+3]);
fclose(f);
signs=[1 1 1 -1 -1 1];
for i=1:6
    uiuj.profile(i,:,:)=0.5*(uiuj.profile(i,:,:)+signs(i)*uiuj.profile(i,:,end:-1:1));
    uiuj.spectra(i,:,:,:)=0.5*(uiuj.spectra(i,:,:,:)+signs(i)*uiuj.spectra(i,:,:,end:-1:1));
    if i<5; uiuj.convs(i,:,:,:)=0.5*(uiuj.convs(i,:,:,:)+signs(i)*uiuj.convs(i,:,:,end:-1:1)); end
end
uiuj.spectra=uiuj.spectra(:,:,:,1:ny/2+2);
uiuj.spectra(:,:,nz+2:end,:)=uiuj.spectra(:,:,nz+2:end,:)+uiuj.spectra(:,:,nz:-1:1,:);
uiuj.spectra=uiuj.spectra(:,:,nz+1:end,:);
uiuj.profile=uiuj.profile(:,:,1:ny/2+2);
uiuj.convs=uiuj.convs(:,:,:,1:ny/2+2);
ttrsp_from_conv=squeeze(sum(uiuj.convs,2));
if premultiply
    uiuj.spectra=uiuj.spectra.*repmat(reshape(kz,[1,1,nz+1,1]),[6,10,1,ny/2+2]);
    ttrsp_from_conv=ttrsp_from_conv.*repmat(reshape(kz,[1,nz+1,1]),[4,1,ny/2+2]);
end


%% Plot various single-point budgets
vars={'(1,1)','(2,2)','(3,3)','(1,2)'};
terms={'var','prod','psdiss','ttrsp','vdiff','pstrain','ptrsp','residual'};
for i=1:4
    figure(i)
    hold on
    title(['budget for ' vars{i} 'component of the Reynolds stress tensor'])
    for j=2:7
        plot(y,squeeze(uiuj.profile(i,j,:)),'linewidth',2)
    end
    plot(y,squeeze(sum(uiuj.profile(i,2:7,:),2)),'linewidth',2)
    legend(terms{2:8})
    xlim([0 1])
end

%% Plot spectral budgets
for i=1:4
    figure(100+i)
    subplot(2,4,1)
    title(['budget for ' vars{i} 'component of the Reynolds stress tensor'])
    hold on
    for k=1:7
       subplot(2,4,k);
       pcolor(lambdaz,y,squeeze(uiuj.spectra(i,k,:,:))'); shading interp
       colorbar;
       title(terms{k});
       set(gca(),'XScale','log','YScale','log');    
    end
    subplot(2,4,8);
    pcolor(lambdaz,y,squeeze(sum(uiuj.spectra(i,2:7,:,:),2))'); shading interp  
    title(terms{8});
    set(gca(),'XScale','log','YScale','log');  
    colorbar();
end

%% Check that convolution give same results
figure(999)
subplot(2,4,1)
for i=1:4
   subplot(2,4,i);
   pcolor(lambdaz,y,squeeze(uiuj.spectra(i,4,:,:))'); shading interp
   colorbar;
   title(['pseudospectral ' terms{4}]);
   set(gca(),'XScale','log','YScale','log');    
   subplot(2,4,4+i);
   pcolor(lambdaz,y,squeeze(ttrsp_from_conv(i,:,:))'); shading interp
   colorbar;
   title(['convolutions ' terms{4}]);
   set(gca(),'XScale','log','YScale','log');    
end

%% Plot triadic interactions
for i=1:5
    figure(200+i)
    hold on
    iy=8;
    if i<5; title(['triadic interactions ' vars{i} ' at height y^+=' num2str(y(iy)*1000)]); end
    iz=find(lambdaz*1000<=(3*(y(iy)*1000)^2),1);
    KZ=(-nz:nz)*beta0; 
    kz0=kz(iz);
    K=1-KZ/kz0; K(K<0)=-K(K<0);
    KZ(KZ<0)=-KZ(KZ<0);
    if i<5
        h=scatter(KZ/kz0,K,[],squeeze(uiuj.convs(i,:,iz,iy))','filled');
    else
        h=scatter(KZ/kz0,K,[],squeeze(sum(uiuj.convs(1:3,:,iz,iy)))','filled');
    end
    cl=get(gca(),'Clim'); cl=max(abs(cl)); set(gca(),'CLim',[-cl cl]);
    colorbar(); colormap redblue;
end