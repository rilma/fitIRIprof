
function N = achapman_feval(h, pars)

    k = 0.5; % alpha term
    
    Nm = pars(1); % maximum electron density (m-3)
    hm = pars(2); % height of maximum electron density (km)
    H = pars(3); % Scale height (km)
    
    z = (h - hm) / H; % Reduced height (km)
    
    N = Nm * exp(k * (1 - z - exp(-z))); % Chapman profile model
        
