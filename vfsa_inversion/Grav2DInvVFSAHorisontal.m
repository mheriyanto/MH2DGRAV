% Author: Mohammad Heriyanto
% Forward Gravity using Talwani Formulation
% ARTICLE: W. Srigutomo, M. Heriyanto, and M. Hilmi Aufa. Gravity Inversion of Talwani Model using Very Fast Simulated Annealing. Journal of Mathematical and Fundamental Sciences, Vol. 51, No. 2, 2019, 177-190. doi: 10.5614/j.math.fund.sci.2019.51.2.7
% URL: journals.itb.ac.id/index.php/jmfs/article/download/8270/4113
% CODE: https://github.com/mheriyanto/MH2DGRAV

clear all; close all; clc;

global kG;
global lpoint;
global z0;
global nx;
global nz;
global M;
global lm;

global si2mg;
global km2m;

global x;
global dGobs;
global rhoobs;
global N;
global m;

global minrho;
global maxrho;
global consVmin;
global consVmax;

%%================ Synthetic Data ==============%%
data = load('syntheticdata.txt');
x = data(:,1);
dGobs = data(:,2);
N = length(x);

%%================== CONSTANT ================%%
kG = 6.673e-11;     % gravity coefisien ( N.m^2/kg^2)
lpoint = 4;
z0 = 0;             % topography

si2mg = 1e5;		    % SI to mGal
km2m = 1e3;			    % km to m

%%================ AUTOMATIC GRID ========================%%
% Square length
dx = x(2)-x(1); 
dz = dx; 
m = dx/2;

% x and z length
panjangx = (x(end)-x(1)) + 2*m; 
panjangz = panjangx/2;

% number of square
nx = panjangx/dx; 
nz = panjangz/dz;

% nz
if mod(nz,1)==0
    nz = nz;
else
    nz = fix(nz);
end

M = nx*nz;
lm = M;

%%================= Block Domain =====================%%

[xm,zm,gz] = DomainBlok(x,dx,dz);

% Change xm and zm domain to j = 1:M
[xf,zf] = ChangeDomainXZ(xm,zm);

%%============== Synthetic Data Plot =================%%
rhoobs = zeros(nz,nx);
for i = 1:nz       % column 
    for j = 1:nx   % raw
        % Horizontal Model
        if i>=4 && i<=7 && j>=7 && j<=14
            rhoobs(i,j) = 1;
        end
        
    end
end

PlotData(x,dGobs,xm,zm,rhoobs)
% print('-dpng','Synthetic Data Plot','-r500');

fprintf('number of column = %i\n',nz);
fprintf('number of raw = %i\n',nx);
fprintf('number of square = %i\n',M);

%%================= Kernell Matrix =====================%%
A = 2*kG.*Kernell(gz);

figure(2)
subplot(1,1,1);
h = pcolor(A);
set(h,'EdgeColor','none');
set(gca,'ydir','reverse');
xlabel(['\bf \fontsize{10}\fontname{Times}The Blocks']);
ylabel(['\bf \fontsize{10}\fontname{Times}Observation point']);
title(['\bf \fontsize{12}\fontname{Times}Matriks A']);

hp = get(subplot(1,1,1),'Position');
hg = colorbar('Position', [hp(1)+hp(3)+0.022  hp(1)+0.022  0.025  hp(2)+hp(3)-0.25]);
title(hg,'Value','fontweight','bold','fontsize',8);
caxis([0 0.3])

% print('-dpng','Plotting Matrix A','-r500');

%%================= Min & Max Model =====================%%
rho = zeros(nz,nx);
minrho = 0;
maxrho = 0.5;

consVmin = 0.35;
consVmax = 1.5;

%%=============== Initial Inversion ======================%%
for i = 1:N
    sumG = 0;
    for j = 1:M
        Gij = A(i,j)*A(i,j);
        sumG = sumG + Gij;
    end
    D(i,i) =  1/sqrt(sumG);
end

lamda = 0.1;
I = eye(N,N);
Vold = A'*D*inv(D*A*A'*D+lamda*I)*D*dGobs;    % Initial Inversion

%%========================== CONSTANT ================================%%
c = 1;
km = 1/lm;
tol = 1e-20;
jj = 1;ii = 1;

% Initial Temperature
qpengali = 100;
T0 = TempAwal(qpengali,A);

T(1) = T0;

itermax = 50;
nv = 20;
tic;

%%================= Intial Forward ============%%
dGcal = A*Vold;
error = norm(dGcal-dGobs);

while (error > tol)
    for (j=1:nv)
      
      % model perturbasi
      %====================== OTHER CONSTRAIN ============================%%
      %%======================== Spesific Random ============================%%
      for i = 1:nz       
          for j = 1:nx   
          rhoit(i,j) = updatemodel(T(jj),Vold(j),maxrho,minrho);
              % Horizontal Model
               if i>=3 && i<=7 && j>=5 && j<=15
                  rhoit(i,j) = updatemodel(T(jj),Vold(j),consVmax,consVmin);
              end
          end
      end

      Vcalit = ChangeDomainRho(rhoit);
      
      %%=============================================================%%
      % next model
      dGcal1 = A*Vold;                            % fwd inital model
      dGcal2 = A*Vcalit;                          % fwd next model

      % misfit value of each model
      E1 = norm(dGcal1-dGobs);
      E2 = norm(dGcal2-dGobs);

      disp('E1 : '); disp(E1);
      disp('E2 : '); disp(E2);

      deltaE(ii) = E2 - E1;                        % misfit difference
      disp('dE : '); disp(deltaE(ii));

      Prob(ii) = exp(-(deltaE(ii))/T(jj));         % probability
      disp('P : '); disp(Prob(ii));

      %%========================= METHOD OF CHOOSING MODEL ======================%%
      if (deltaE(ii) <= 0); 
          Vold = Vcalit;
          dGCal = dGcal2; 
          error = E2;
      elseif (deltaE(ii) > 0);
          bil_rand = rand(1,1);
          if (Prob(ii)> bil_rand);
              disp('cek');
              Vold = Vcalit;
              dGCal = dGcal2; 
              error = E2;
          else
              Vold = Vold;
              dGCal = dGcal1; 
              error = E1;
          end
      end

      %%======== Iteration Print =========%%
      fprintf('==================\n');
      fprintf('Error = %f\n',error);
      fprintf('Iteration = %i\n',ii);
      % Error RMS
      RMSE(ii) = sqrt(error)/N;

      %============================ ITERATION PLOT ===============================%
      rho = DeChangeRho(Vold);            % Change domain (Vold to rho)
      
      PlotFig(rho,dGcal,xm,zm,ii,error);
      % if(ii==1 || (ii >= 1 && mod(ii,100)==0))
      %      saveas(gcf,['Inversion VFSA Iteration-',num2str(ii),'.png']);
      % end

      iterplot(ii) = ii;
      misfit_plot(ii) = error;
      ii = ii+1;
    end

    %===================== STOPPING CRITERION =======================%
    iter(jj) = jj;
    if(iter(jj) > itermax)
        break
    end

    s = 0;
    for j = 1:M
        if Vold(j) <= Vold(j)
           s = s + 1;  
        end
    end
    fprintf('misfit = %f || s = %i \n',error,s);
    %     if s == M
    %         break;
    %     end

    %============================ TEMPERATURE decrease =======================%
    T(jj+1) = T(jj)*exp(-c*(jj^(km)));
    jj = jj+1;
end
toc


figure();
subplot(2,1,1);
plot(x,dGobs,'*','color','r','MarkerSize',5);
hold on
plot(x,dGcal,'-','color','b','LineWidth',2);
xlabel(['\bf \fontsize{10}\fontname{Times}Distance (m)']);
ylabel(['\bf \fontsize{10}\fontname{Times}\Delta g (mGal)']);
leg = legend('Obs','Cal'); set(leg,'Location','South','fontsize',8);

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

ylabel(['\bf \fontsize{10}\fontname{Times}Depth (m)']);
set(gca,'ydir','reverse');
xlim([min(x) max(x)]);
axis tight;

hp = get(subplot(2,1,2),'Position');
% bla,atas,
hcb = colorbar('Position', [hp(1)+hp(3)+0.022  hp(1)+0.022  0.025  hp(2)+hp(3)-0.25]);
title(hcb,'g/cm^3','fontweight','bold','fontsize',10);
caxis([0 1]);

print('-dpng','Plot Final VFSA','-r500');

%======================= SAVING PARAMETERS ==========================%
figure(4)
subplot(2,2,[1 2])
semilogx(iterplot,misfit_plot,'b','LineWidth',2);
xlabel('Iteration','fontweight','bold','fontsize',11);
ylabel('Misfit','fontweight','bold','fontsize',11);
axis tight
grid on;

subplot(2,2,[3 4])
plot(iter,T,'b','LineWidth',2);
xlabel('Iteration','fontweight','bold','fontsize',11);
ylabel('Temperature','fontweight','bold','fontsize',11);
axis tight
grid on;

print('-dpng','Plot Parameter VFSA 1','-r500');

figure(5)
subplot(2,2,[1 2])
semilogx(iterplot,RMSE,'b','LineWidth',2);
xlabel('Iteration','fontweight','bold','fontsize',11);
ylabel('RMSE','fontweight','bold','fontsize',11);
axis tight
grid on;

subplot(2,2,[3 4])
plot(iter,T,'b','LineWidth',2);
xlabel('Iteration','fontweight','bold','fontsize',11);
ylabel('Temperature','fontweight','bold','fontsize',11);
axis tight
grid on;

print('-dpng','Plot Parameter VFSA 2','-r500');

%============================ SAVE WORKSPACE OF RESULTS ===============================%

save Grav2DInvVFSA.mat