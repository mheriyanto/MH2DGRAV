function [xf,zf] = ChangeDomainXZ(xm,zm)
%CHANGEDOMAINXM Summary of this function goes here
%   Detailed explanation goes here

  global nz;
  global nx;
  global m;

  k = 1;
  for i = 1:nz       % coloumn
      for j = 1:nx   % raw
          xf(k,1) = xm(i,j) + m;
          zf(k,1) = zm(i,j) + m;
          k = k + 1;
      end
  end

end