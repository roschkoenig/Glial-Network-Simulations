function dxdt = glia_cortsheet(t,x,T)

% Initialise 
%--------------------------------------------------------------------------
dxdt    = zeros(length(x),1);
alli    = 1:length(x);
thisi   = [1:fix(alli(end)/3)]*3-2;
E      = thisi;     % excitatory populations
I      = thisi+1;   % inhibitory populations
G      = thisi+2;   % glia cell populations

% Connectivity parameters
%==========================================================================
% Excitatory to excitatory connectivity: extrinsic and intrinsic
%--------------------------------------------------------------------------
if ~isfield(T, 'c_EE') 
    c_EE    = randn(length(E)) * 0;
    for c = 1:length(c_EE)
        c_EE(c,c) = 24;
    end
else
    c_EE = T.c_EE;
end

% Excitatory to inhibitory connectivity: intrinsic
%--------------------------------------------------------------------------
if ~isfield(T, 'c_EI')
    c_EI    = zeros(length(E));
    for c = 1:length(c_EI)
        c_EI(c,c) = 40;
    end
else
    c_EI = T.c_EI;
end

% Inhibitory to excitatory connectivity: intrinsic
%--------------------------------------------------------------------------
if ~isfield(T, 'c_IE')
    c_IE    = zeros(length(E));
    for c = 1:length(c_IE)
        c_IE(c,c) = -20;
    end
else
    c_IE = T.c_IE;
end

% Inhibitory to inhibitory connectivity: intrinsic
%--------------------------------------------------------------------------
if ~isfield(T, 'c_II')
    c_II    = zeros(length(E));
    for c = 1:length(c_II)
        c_II(c,c) = 0;
    end
else
    c_II = T.c_II;
end

% Excitatory to glia connectivity: intrinsic
%--------------------------------------------------------------------------
if ~isfield(T, 'c_II')
    c_II    = zeros(length(E));
    for c = 1:length(c_II)
        c_II(c,c) = 0;
    end
else
    c_II = T.c_II;
end

% Glia to excitatory connectivity: intrinsic
%--------------------------------------------------------------------------
if ~isfield(T, 'c_GE')
    c_GE    = zeros(length(E));
    for c = 1:length(c_GE)
        c_GE(c,c) = 0;
    end
else
    c_GE = T.c_GE;
end

% Excitatory to glia connectivity: intrinsic
%--------------------------------------------------------------------------
if ~isfield(T, 'c_EG')
    c_EG    = zeros(length(E));
    for c = 1:length(c_EG)
        c_EG(c,c) = 0;
    end
else
    c_EG = T.c_EG;
end



% Other model parameters
%==========================================================================
% Baseline activation
%--------------------------------------------------------------------------
if ~isfield(T, 'P'), P = ones(length(E),1) * (1.5);  	else, 	P = T.P;    end;  
if ~isfield(T, 'Q'), Q = ones(length(E),1) * (-2);      else,   Q = T.Q;    end;

% Time constants
%--------------------------------------------------------------------------
tau(1,1)    = .013;
tau(1,2)    = .013;
tau(1,3)    = .13;

% Contribution of random noise
%--------------------------------------------------------------------------
if ~isfield(T, 'A_S'),  A_S = 0;    else,   A_S = T.A_S;    end

%==========================================================================
% ODE system
%==========================================================================
dxdt(E) = ( - x(E) + sig(c_EE * x(E) + c_IE * x(I) + c_GE * x(G) + P + A_S * randn(length(E),1)) ) / tau(:,1);
dxdt(I) = ( - x(I) + sig(c_EI * x(E) + c_II * x(I) + Q) ) / tau(:,2);
dxdt(G) = ( - x(G) + c_EG * x(E)) / tau(:,3);
 
end


% Sigmoid function
%--------------------------------------------------------------------------

function z = sig(u)

    alpha   = 1;
    theta   = 4;
    
    z = 1 ./ (1 +exp( - alpha * (u - theta )));
end     
