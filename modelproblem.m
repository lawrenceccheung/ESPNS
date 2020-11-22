% modelproblem.m
% ----
%
% Create the polynomial system for the model problem:
%     u'' + u^2 - f(x) = 0
% for 0 <= x < 2*pi, with periodic boundary conditions:
%     u(0) = u(2*pi)   and   f(0) = f(2*pi)
%
clear all
addpath("PNLA_MATLAB_OCTAVE");
addpath("SuiteSparse/SPQR/MATLAB")

% -- Translate k index to matlab index
ktoi = @(k,kmin) (k-kmin+1);
% --

% Set the size of the system
kmin = -3;
kmax = +3;
%stol=0.5E-14;  % |k| <= 2
stol=2.5E-14;  % |k| <= 3


% Set the f vector
f    = zeros(1,length(kmin:kmax));
f(ktoi(1, kmin)) = 1;

% A vector with all k integers
allk = kmin:kmax;

% The polynomial system
polysys = cell(length(allk), 2);

% Construct each k equation
for k=allk
    eqncoeffs = [];
    kpowers   = [];
    
    % Add the -k^2 term
    if (k ~= 0)
        eqncoeffs(end+1) = -k*k;
        upower = 0*allk;
        upower(ktoi(k, kmin)) = 1;
        kpowers(end+1,:)     = upower;
    end
    
    % Add the Cpqk term
    for p=allk
        q = k-p;  % Cpqk
        if ((kmin <= q) && (q <= kmax))
            upower = 0*allk;
            if (p<q)
                eqncoeffs(end+1) = 2;
                upower(ktoi(p, kmin)) = 1;
                upower(ktoi(q, kmin)) = 1;
                kpowers(end+1,:) = upower;
            elseif (p == q)
                eqncoeffs(end+1) = 1;
                upower(ktoi(p, kmin)) = 2;
                kpowers(end+1,:) = upower;
            end
        end
    end
    
    % Add the f constant term
    if abs(f(ktoi(k,kmin)) > 1.0E-16)
        eqncoeffs(end+1) = f(ktoi(k,kmin));
        upower = 0*allk;        
        kpowers(end+1,:) = upower;    
    end
    
    % Add things to the polysys cell matrix
    polysys{ktoi(k,kmin), 1} = eqncoeffs;
    polysys{ktoi(k,kmin), 2} = kpowers;
end

% Now solve the polynomial system
% ------
tic;

%[root, d, c, ns, check, cr, digits] = sparf2(polysys);
[root, d, c, ns, check, cr, digits] = qdsparf2(polysys, stol);

elapsedtime = toc;
calcerror(root, kmin, kmax);
fprintf('Solve time: %f sec\n',elapsedtime)
