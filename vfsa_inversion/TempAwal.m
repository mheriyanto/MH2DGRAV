function T = TempAwal(qpengali,A)
%TEMPAWAL Summary of this function goes here
%   Detailed explanation goes here

  global M;
  global minrho;
  global maxrho;
  global dGobs;

  for aa = 1:1:30
    for dd = 1:1:M
        Vtem(dd,1) = minrho + rand(1,1)*(maxrho-minrho);
    end
    dGcal = A*Vtem;
    misfit_tem = norm((dGcal-dGobs)/dGobs);
    misfit_tem_accum(aa) = misfit_tem;
  end

  stdTemp = std(misfit_tem_accum);
  T = stdTemp*qpengali;                % Initial temperature

end