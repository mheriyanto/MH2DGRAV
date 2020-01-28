% Mohammad Heriyanto
% Forward Gravity using Talwani Formulation
% W. Srigutomo, M. Heriyanto, and M. Hilmi Aufa. Gravity Inversion of Talwani Model using Very Fast Simulated Annealing. Journal of Mathematical and Fundamental Sciences, Vol. 51, No. 2, 2019, 177-190. doi: 10.5614/j.math.fund.sci.2019.51.2.7
% URL: https://github.com/mheriyanto/MH2DGRAV

clear all; close all; clc;

global kG;
global lpoint;
global z0;
global nx;
global nz;
global M;

global si2mg;
global km2m;

global N;
global m;

%=== Constants ===%
kG = 6.673e-11;     % Gravity const, N.m^2/kg^2
lpoint = 4;
z0 = 0;             % topography

si2mg = 1e5;		% SI to mGal
km2m = 1e3;			% km to m

%=== GRID ===%
% space dimension
panjangx = 1000; panjangz = 500;

% box length
dx = 50; 
dz = dx; 
m = dx/2;

% number of box
nx = panjangx/dx; 
nz = panjangz/dz;
M = nx*nz;

fprintf('coloum number = %i\n',nz);
fprintf('row number = %i\n',nx);
fprintf('box number = %i\n',M);

%==== location of measurement ===%
for k = 1:nx
    x(k) = dx*k - m; 
end
N = length(x);

%=== initial model (density) ===%
rho = zeros(nz,nx);
for i = 1:nz       % Kolom
    for j = 1:nx   % Baris
        
        % horizontal model
%          if i>=4 && i<=7 && j>=7 && j<=14
%              rho(i,j) = 1;
%          end
        
        % vertical model
        if i>=2 && i<=8 && j>=10 && j<=11
            rho(i,j) = 1;
        end
        
    end
end

% arrange Rho to V matrix
V = ChangeDomainRho(rho);

%=== Block Domain ===%
[xm,zm,gz] = DomainBlok(x,dx,dz);

%=== Kernell matrix ===%
A = 2*kG.*Kernell(gz);

figure(2)
subplot(1,1,1)
h = pcolor(A);
set(h,'EdgeColor','none');
set(gca,'ydir','reverse');
xlabel(['\bf \fontsize{14}\fontname{Times}The Blocks']);
ylabel(['\bf \fontsize{14}\fontname{Times}Observation point']);

hp = get(subplot(1,1,1),'Position');
% bla,atas,
hg = colorbar('Position', [hp(1)+hp(3)+0.022  hp(1)+0.022  0.025  hp(2)+hp(3)-0.25]);
title(hg,'Value','fontweight','bold','fontsize',8);
caxis([0 0.3])

print('-dpng','Matrix A Plot','-r500');

%=== Forward Calculation ===%
dGcal = A*V;

%=== Add Noise to Data ===%
for k=1:N
    noise(k) = 0.1*dGcal(k);              %noise 10%
    dGcal_noise(k) = dGcal(k) + noise(k)*(rand()-0.5);
end

%=== Inversion Parameters ===%
figure(3);
subplot(2,1,1);
plot(x,dGcal_noise,'*','color','r','MarkerSize',5);
xlabel(['\bf \fontsize{12}\fontname{Times}Distance (m)']);
ylabel(['\bf \fontsize{12}\fontname{Times}\Delta g (mGal)']);
xlim([min(x)-m max(x)+m]);

subplot(2,1,2);
hold on;

for i=1:nz
    for j= 1:nx
        xleft(i,j) = xm(i,j);
        xright(i,j) = xm(i,j+1);
        xfill = [xleft(i,j) xright(i,j) xright(i,j) xleft(i,j)];
        zfill = [zm(i,j) zm(i,j+1) zm(i+1,j+1) zm(i+1,j)];
        color = rho(i,j);
        patch(xfill,zfill,color,'EdgeColor','none');
    end
end

ylabel(['\bf \fontsize{12}\fontname{Times}Depth (m)']);
set(gca,'ydir','reverse');
xlim([min(x) max(x)]);
axis tight;

hp = get(subplot(2,1,2),'Position');
hcb = colorbar('Position', [hp(1)+hp(3)+0.022  hp(1)+0.022  0.025  hp(2)+hp(3)-0.25]);
title(hcb,'g/cm^3','fontweight','bold','fontsize',10);

print('-dpng','Plot forward of vertical model','-r500');
simpan = [x' dGcal_noise'];
save('DataVerticalBox.txt','simpan','-ascii');
