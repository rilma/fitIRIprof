
function [m, nloop] = lm_method(funcname, jacname, m0, max_nloop, ydata, sigma)
%
% Applies The Levenberg-Marquardt method to solve a nonlinear problem in a
% least-squares sense
%

    % Modification factor of initial guess and 2-norm of the square errors
    % normalized to their respective variance
    factor = 2.0;
    
    % Termination tolerance
    tolerance = 1e-12;
    
    % Positive parameter to be adjusted to ensure convergence
    lambda = 0.0001;    
    
    % Modification factors of the "lambda" positive parameter
    lfactorUP = 2.5; lfactorDOWN = 0.5;
    
    % The "lambda" postive parameter can be equal to values within the
    % interval pointed out by these constants:
    minLambda = 1e-12; maxLambda = 1e16;
    
    % Selecting the initial guess of model parameters
    m = m0;
    
    % Implementation of eq. 9.13
    %
    
    % errors normalized by their respective standard deviation   
    fvFunc = (feval(funcname, m) - ydata) ./ sigma;
    
    % Sum of the squared errors normalized by their respective standard
    % deviation
    fm = norm(fvFunc, 2) ^ 2;
    
    %
    
    % Assigning the modified inital guess and 2-norm to a temporal
    % variables utilized in the iterative process below
    mTMP = m0 * factor; fmTMP = fm * factor;
    
    % Reshaping the "sigma" vector to be utilized with the Jacobian values
    sigmaNd = reshape(repmat(sigma, numel(m0), 1), numel(sigma), numel(m0));    
    
    
    % Iterations
    
    nloop = 0;
    
    while (nloop <= max_nloop)
        
        fvJac = feval(jacname, m);        
        
        J = fvJac ./ sigmaNd;
        
        fvFunc = (feval(funcname, m) - ydata) ./ sigma;
    
        % RHS of eq. 9.30 (Levenberg-Marquardt method)
        jm = -J' * fvFunc;
                
        jfTMP = [-fvFunc; zeros(length(m), 1)];
        
        dm = [J; sqrt(lambda) * eye(length(m))] \ jfTMP;
        
        % Termnination criteria
        if ~(((norm(jm, 2) < sqrt(tolerance) * (1 + abs(fm))) & ...                
                (norm(mTMP - m, 2) < sqrt(tolerance) * (1 + norm(m, 2)))) & ...
                (abs(fmTMP - fm) < tolerance * (1 + abs(fm))))
            
            % Upgrading value of eq. 9.13        
            %
            fvFuncS = (feval(funcname, m + dm) - ydata) ./ sigma;
        
            fmUP = norm(fvFuncS, 2) ^ 2;
        
            %

            if (fmUP < fm)
            
                mTMP = m; m = m + dm;
                fmTMP = fm; fm = fmUP;
                lambda = lambda * lfactorDOWN;
            
                if lambda < minLambda, lambda = minLambda; end;
            
            else
            
                lambda = lambda * lfactorUP;
            
                if lambda > maxLambda, lambda = maxLambda; end;
            
            end
        
            nloop = nloop + 1;
        
        else
            
            % it quit the loop and function when termination criteria is
            % not satisfied
            return;
            
        end
        
    end
