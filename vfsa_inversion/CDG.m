function Ab = CDG(gz)
%CHANGEDOMAIN Summary of this function goes here
%   Detailed explanation goes here

  global nz;
  global nx;

  % Change Domain
  k = 1;
  for i = 1:nz
     for j = 1:nx
         Ab(k) = gz(1,i,j);
         k = k + 1;
     end
  end
end