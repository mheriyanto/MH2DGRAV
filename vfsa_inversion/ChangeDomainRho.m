function V = ChangeDomainRho(rho)
%CHANGEDOMAINRHO Summary of this function goes here
%   Detailed explanation goes here

  global nz;
  global nx;
  k = 1;
  for i = 1:nz       % coloumn
      for j = 1:nx   % raw
          V(k,1) = rho(i,j);
          k = k + 1;
      end
  end

end