function [xm,zm,gz] = DomainBlok(x,dx,dz)
%DOMAINBLOK Summary of this function goes here
%   Detailed explanation goes here

  global nx;
  global nz;
  global N;
  global z0;
  global lpoint;

  % Filling value of each elements
  for i = 1:nz+1
      for j = 1:nx+1
          xm(i,j) = (j-1)*dx; 
          zm(i,j) = (i-1)*dz;
      end
  end

  for n = 1:N
      for i = 1:nz
          for j = 1:nx
              % 4 points of xm and zm
              xmm = [xm(i,j);xm(i,j+1);xm(i+1,j+1);xm(i+1,j)];
              zmm = [zm(i,j);zm(i,j+1);zm(i+1,j+1);zm(i+1,j)];
          
              % One Block ( 4 points & 1 rho)
              gz(n,i,j) = Talwani(x(n),z0,xmm,zmm,lpoint); 
          end
      end
  end

  end