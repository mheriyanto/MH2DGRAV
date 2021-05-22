function A = Kernell(gz)
%KERNELL Summary of this function goes here
%   Detailed explanation goes here

  global N;
  global M;

  for n = 1:N
      Ab = CDG(gz(n,:,:));
      Ac{n} = Ab;
  end

  % Anomaly in 1 point of measurement           
  for k = 1:N
     for l = 1:M
         A(k,l) = Ac{k}(l);
     end
  end

end