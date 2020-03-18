function [dns]=init_dns(filename,walls)

f = fopen(filename,'r');
l=fgetl(f); l=str2num(l(1:find(l=='!')-1));
dns.nx=l(1); dns.ny=l(2); dns.nz=l(3);
l=fgetl(f); l=str2num(l(1:find(l=='!')-1));
dns.alfa0=l(1); dns.beta0=l(2);
l=fgetl(f); l=str2num(l(1:find(l=='!')-1));
dns.Re=l(1); dns.nu=1/dns.Re;
l=fgetl(f); l=str2num(l(1:find(l=='!')-1));
dns.a=l(1); dns.ymin=l(2); dns.ymax=l(3);
l=fgetl(f); l=str2num(l(1:find(l=='!')-1));
dns.meanpx=l(1); dns.meanpz=l(2);
l=fgetl(f); l=str2num(l(1:find(l=='!')-1));
dns.meanflowx=l(1); dns.meanflowz=l(2);
l=fgetl(f); l=str2num(l(1:find(l=='!')-1));
dns.deltat=l(1); dns.cflmax=l(2); dns.time=l(3);

dns.nxd = floor(3*dns.nx/2); while ~fftfit(dns.nxd); dns.nxd=dns.nxd+1; end
dns.nzd = 3*dns.nz         ; while ~fftfit(dns.nzd); dns.nzd=dns.nzd+1; end

if walls==2
    dns.y = 0.5*(1 + tanh( dns.a * ( 2 * ( (-1:dns.ny+1) )/dns.ny - 1) )/tanh( dns.a ))*(dns.ymax-dns.ymin) + dns.ymin;
elseif walls==1
    dns.y = (1 + tanh( dns.a * ( ( (-1:dns.ny+1) )/dns.ny - 1) )/tanh( dns.a ))*(dns.ymax-dns.ymin) + dns.ymin;
end
dns.x = (0:2*dns.nxd-1)*(2*pi/dns.alfa0)/(2*dns.nxd);
dns.z = (0:dns.nzd-1)*(2*pi/dns.beta0)/dns.nzd;
dns.size = [dns.ny+3,2*dns.nz+1,dns.nx+1,3];
dns.sized = [dns.ny+3,dns.nzd,2*dns.nxd,3];

end 