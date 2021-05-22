function mcalit = updatemodel(T,mcal,maxm,minm)
%UPDATEMODEL Summary of this function goes here
%   Detailed explanation goes here

    u = rand(1,1);
    yu = sign(u-0.5)*T*(((1+(1/T))^(abs(2*u-1)))-1);
    mcalit = mcal + yu*(maxm-minm);
    % kondisi rhomin <= rho <= rhomax
    while (mcalit < minm || mcalit > maxm) 
       u = rand(1,1);
       yu = sign(u-0.5)*T*(((1+(1/T))^(abs(2*u-1)))-1);
       mcalit = mcal + yu*(maxm-minm);
    end        
end

