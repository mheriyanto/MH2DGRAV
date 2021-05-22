function rho = DeChangeRho(V)
%DECHANGERHO Summary of this function goes here
%   Detailed explanation goes here

  global nz;
  global nx;

  k = 1;
  for i = 1:nz       % coloum
      for j = 1:nx   % raw
          rho(i,j) = V(k,1);
          k = k + 1;
      end
  end

end