function g = Talwani(x0,z0,xcorn,zcorn,ncorn)

global si2mg;
global km2m;

% J.A. Olive, MIT/WHOI Joint Program
% Potential Theory in Gravity and Magnetic Applications, Blakely, 1995, p. 378
% all distance parameters in units of km, rho in units of kg/m^3
% Output = vertical attraction of gravity g, in mgal

% run
sumG = 0;
for n = 1:ncorn		% number of angle
    if n == ncorn
        n2 = 1;
    else
        n2 = n+1;
    end
    
    % calculate one side of polygon, (x1,z1) and (x2,z2)
    x1 = xcorn(n)-x0;
    z1 = zcorn(n)-z0;
    x2 = xcorn(n2)-x0;
    z2 = zcorn(n2)-z0;
	
    % calculate r1 and r2
    r1sq = x1*x1 + z1*z1;
    r2sq = x2*x2 + z2*z2;
    
    % condition r1 or r2 = 0 
    if r1sq == 0
        break;
        disp('Field point on corner');
    elseif r2sq == 0
        break;
        disp('Field point on corner');
    end
    
    % condition z2-z1 = 0, so the difference must be made of value
    denom = z2 - z1;
    if denom == 0
        denom = 1e-6;
    end
   
    alpha = (x2-x1)/denom; 	% alfa = (x2-x1)/(z2-z1)
    beta = (x1*z2-x2*z1)/denom; % beta = x1-alfa*z1
	
    % first part
    factor = beta/(1+alpha*alpha);
	
    % second part
    term1 = 0.5*(log(r2sq)-log(r1sq));
    term2 = atan2(z2,x2)-atan2(z1,x1);
	
    % sum of all
    sumG = sumG + factor*(term1-alpha*term2);
end

% satu bentuk loop tertutup
g = sumG*si2mg*km2m;

end
