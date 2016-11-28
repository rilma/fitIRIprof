
function achjc = achapman_jac(pars)

    global hvalues;
    
    k = 0.5; % alpha term
    
    Nm = pars(1); % maximum electron density (m-3)
    hm = pars(2); % height of maximum electron density (km)
    H = pars(3); % Scale height (km)
    
    z = (hvalues - hm) / H; % Reduced height (km)
    
    %
    %
    
    % Initiali...
    achjc = repmat(NaN, numel(hvalues), 3);
    
    % 1st derivative respect to Nm
    achjc(:, 1) = achapman_funct(pars) / Nm;
    
    % 1st derivative respect to hm
    achjc(:, 2) = achapman_funct(pars) .* ((k / H) * (1 - exp(-z)));
    
    % 1st derivative respect to H
    achjc(:, 3) = achapman_funct(pars) .* ((k * z / H) .* (1 - exp(-z)));
    