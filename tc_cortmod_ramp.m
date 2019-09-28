function B = tc_cortmod_ramp

% Housekeeping
%--------------------------------------------------------------------------
tc_housekeeping;
rng(45)

% Define parameters of interest manually
%==========================================================================
% All nodes
%--------------------------------------------------------------------------
Sarray      = 25; 
Ns          = Sarray^2;    
A_S         = 0;    % scales contribution of random input to activity
connected   = 1;    % builds a connectivity matrix between nodes; 0 = unconnected

Plim       =  [1 1.5];    % Maximum excitability of the system
Ptype      = 'randwalk';

% Standard nodes
%--------------------------------------------------------------------------
E  	= tc_modelparams('hopf');

EE  = E.EE;             EI  = E.EI;
II  = E.II;             IE  = E.IE;
P   = E.P;              Q   = E.Q; 

% Epileptic node
%--------------------------------------------------------------------------
epi     = 1;         	% Switches special epileptic node on or of
e       = 1;            % Index of the abnormal node: fix(Ns / 2);

eEE  = E.EE;            eEI  = 30;
eII  = E.II;            eIE  = E.IE;
eP   = E.P;           	eQ   = E.Q; 

% ODE options
%--------------------------------------------------------------------------
options     = odeset('MaxStep',0.5);
t_range     = 1:.5:10000;
x_ini       = rand(Ns * 2,1) / 10;  % 0 -> near lower fix point; 

% Generate time course for P
%--------------------------------------------------------------------------
switch Ptype 
    case 'ramp'
        tquart  = floor(length(t_range)/4);
        tslope  = range(Plim)/tquart; 
        for t = 1:tquart,                       P1(t)   = Plim(1);  end
        for t = [1:tquart] + tquart,            P1(t)   = tslope * (t-tquart) + Plim(1);  end
        for t = [1:tquart] + 2*tquart,          P1(t)   = Plim(2); end
        for t = 3*tquart+1 : length(t_range),   P1(t) = Plim(2) - tslope * (t - 3*tquart); end

        T.P         = repmat(P1,Ns,1);
        
    case 'randwalk'
        p   = randn(Ns,length(t_range)+500);
        np  = zeros(size(p,1), length(t_range));
        for pp = 1:size(p,1)
            tp      = p(pp,:);
            tp = smooth(tp,500);  
            tp = tp(250:end-251);
            tp = tp / range(tp) * range(Plim) + mean(Plim);
            np(pp,:) = tp;
        end
        
        T.P = np;
        
end
T.tlist     = t_range;


% Run simulations across different parameter values
%==========================================================================
clear B

% Convert parameter values for use in model function
%==================================================================
% Excitatory connections within and between sources
%------------------------------------------------------------------
if connected    

    % Creates connectivity matrix A
    %--------------------------------------------------------------
    i           = 0; 

    for nr = 1:Sarray
    for nc = 1:Sarray
        i   = i + 1;
        xloc(i,:) = [nr,nc];

    end
    end
    D = dist(xloc');

    mu          = 0;
    sg          = 2;
    rdfact      = 0.2;      % Scales additional random connectivity
    A           = tc_dist2cnx(D, mu, sg, rdfact, Ns);
    A           = [ Ns * EE / sum(sum(A)) ] * A;

    T.c_EE      = A;

else 
    c_EE    = zeros(Ns); 
    for c = 1:length(c_EE),     c_EE(c,c) = EE;     end
    T.c_EE  = c_EE;
end

% Further within source connectivity
%------------------------------------------------------------------
c_EI    = zeros(Ns); 
for c = 1:length(c_EI),     c_EI(c,c) = EI;     end
T.c_EI  = c_EI;

c_II    = zeros(Ns); 
for c = 1:length(c_II),     c_II(c,c) = II;     end
T.c_II  = c_II;

c_IE    = zeros(Ns); 
for c = 1:length(c_IE),     c_IE(c,c) = IE;     end
T.c_IE  = c_IE;

% Other parameters
%------------------------------------------------------------------
T.A_S       = A_S;
T.Q         = ones(Ns,1) * Q;

% Special 'Epileptic' Node
%------------------------------------------------------------------
if epi
    if connected,   meanC       = sum(sum(A)) / (2 * Ns);
                    T.c_EE(e,e) = meanC * eEE / EE;
%                     T.c_EE(:,e) = 1.5 * T.c_EE(:,e);    % Rescale connectivity from epileptic node
%                     T.c_EE(e,:) = 0.5 * T.c_EE(e,:);    % Rescale connectivity to epileptic node

    else,           T.c_EE(e,e)     = eEE;  
    end

    T.c_EI(e,e)     = eEI;      
    T.c_II(e,e)     = eII;
    T.c_IE(e,e)     = eIE;
end


% Run actual simulations
%------------------------------------------------------------------
[t,x]       = ode45(@tc_cortsheet, t_range, x_ini, options, T);

% Store in output structure
%------------------------------------------------------------------
B.x      = x;
B.t      = t;
B.T      = T;
B.P      = T.P;

subplot(2,1,1)
    plot(B.t, B.x)
    
subplot(2,1,2) 
    plot(B.t, B.P);
    
end

%==========================================================================
% Model definition
%==========================================================================

function dxdt = tc_cortsheet(t,x,T)

% Initialise 
%--------------------------------------------------------------------------
dxdt    = zeros(length(x),1);
alli    = 1:length(x);
E      = alli(mod(alli,2) == 1);   % odd rows      = excitatory populations
I      = alli(mod(alli,2) == 0);   % even rows     = inhibitory populations

% Connectivity parameters
%==========================================================================
% Connectivity: extrinsic and intrinsic
%--------------------------------------------------------------------------
c_EE = T.c_EE;
c_EI = T.c_EI;
c_IE = T.c_IE;
c_II = T.c_II;

% Other model parameters
%==========================================================================
% Time varying P parameter
%--------------------------------------------------------------------------
tlist       = T.tlist;
[c index]   = min(abs(tlist-t));
P           = T.P(:,index);

% Background inhibition
%--------------------------------------------------------------------------
if ~isfield(T, 'Q'), Q = ones(length(E),1) * (-2);      else,   Q = T.Q;    end;

% Time constants
%--------------------------------------------------------------------------
tau(1,1)    = .013;
tau(1,2)    = .013;

% Contribution of random noise
%--------------------------------------------------------------------------
if ~isfield(T, 'A_S'),  A_S = 0;    else,   A_S = T.A_S;    end

%==========================================================================
% ODE system
%==========================================================================
dxdt(E) = ( - x(E) + sig(c_EE * x(E) + c_IE * x(I) + P + A_S * randn(length(E),1)) ) / tau(:,1);
dxdt(I) = ( - x(I) + sig(c_EI * x(E) + c_II * x(I) + Q) ) / tau(:,2);

 
end


% Sigmoid function
%--------------------------------------------------------------------------

function z = sig(u)

    alpha   = 1;
    theta   = 4;
    
    z = 1 ./ (1 +exp( - alpha * (u - theta )));
end     

% B = tc_cortmod_ramp;
% 
% %%
% dat = ft_preproc_bandpassfilter(B.x', 1000, [1 200]);
% E   = 1:2:size(dat,1);
% 
% for e = 1:length(E)
%     ee = E(e);
%     plot(dat(ee,:)+e/4); hold on
% end