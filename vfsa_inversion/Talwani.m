function g = Talwani(x0,z0,xcorn,zcorn,ncorn)

  global si2mg;
  global km2m;

  % J.A. Olive, MIT/WHOI Joint Program
  % Potential Theory in Gravity and Magnetic Applications, Blakely, 1995, p. 378
  % all distance parameters in units of km, rho in units of kg/m^3
  % Output = vertical attraction of gravity g, in mgal
  % URL: https://jaolive.weebly.com/codes.html

  % run
  sumG = 0;
  for n=1:ncorn		        % banyaknya sudut
    if n == ncorn
        n2 = 1;
    else
        n2 = n+1;
    end
    
    % menghitung 1 sisi polygon, (x1,z1) dan (x2,z2)
    x1 = xcorn(n)-x0;
    z1 = zcorn(n)-z0;
    x2 = xcorn(n2)-x0;
    z2 = zcorn(n2)-z0;
  
    % menghitung r1 dan r2
    r1sq = x1*x1 + z1*z1;
    r2sq = x2*x2 + z2*z2;
    
    % kondisi ketika r1 atau r2 = 0 
    if r1sq == 0
        break;
        disp('Field point on corner');
    elseif r2sq == 0
        break;
        disp('Field point on corner');
    end
    
    % kondisi dimana z2-z1 = 0, maka harus dibuat selisihnya ada nilainya
    denom = z2 - z1;
    if denom == 0
        denom = 1e-6;
    end
   
    alpha = (x2-x1)/denom; 		    % alfa = (x2-x1)/(z2-z1)
    beta = (x1*z2-x2*z1)/denom;   % beta = x1-alfa*z1
  
    % suku pertama
    factor = beta/(1+alpha*alpha);
  
    % suku kedua
    term1 = 0.5*(log(r2sq)-log(r1sq));
    term2 = atan2(z2,x2)-atan2(z1,x1);
  
    % total semua
    sumG = sumG + factor*(term1-alpha*term2);
  end

  % satu bentuk loop tertutup
  g = sumG*si2mg*km2m;

end