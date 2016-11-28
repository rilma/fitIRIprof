
function [N, J] = iri_f2bottom_functeval(m, h)

    global Nm hm;

    B0 = m(1); B1 = m(2);
    
    x = (hm - h) / B0;
    
    N = Nm * exp(-x .^ B1) ./ cosh(x);

    if nargout > 1
        
    % Initiali...
        J = repmat(NaN, numel(h), 2);
    
    % 1st derivative respect to B0
        J(:, 1) = N .* (B1 / B0) .* x .^ B1 .* ...
            (1 + tanh(x) / B1);
  
    % 2nd derivative respect to B1
        J(:, 2) = -N .* x .^ B1 .* log(x);
    
    end
