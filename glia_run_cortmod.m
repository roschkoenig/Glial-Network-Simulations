% Housekeeping
%--------------------------------------------------------------------------
clear F
fs          = filesep;
F.base      = '/Users/roschkoenig/Dropbox/Research/1909 CB Dispatch Sims/Network Sims';
F.scripts   = [F.base fs 'Scripts'];
addpath(genpath(F.scripts));

rng(45)


% Define parameters of interest manually
%==========================================================================
% All nodes
%--------------------------------------------------------------------------
Ns          = 10;
A_S         = 0.1;    % scales contribution of random input to activity
connected   = 1;    % builds a connectivity matrix between nodes; 0 = unconnected
dir         = 2;    % 1 = Forward only; 2 = Forward and Backward
steps       = 150;    % How many simulation steps
sync        = 1;
Pvary       = 'IE';         
Prange      = linspace(-80, 0, steps); 


clear C
for cond = 1:2
% Standard nodes
%--------------------------------------------------------------------------
EE  = 23;               EI  = 35;
II  = 0;                IE  = 15;
P   = .5;               Q   = E.Q;    
if cond == 1,           EG  = 20;    GE  = 1;   % [0,0] vs [20,1] = no glia vs glia
else,                   EG  = 0;     GE  = 0;   end % Switch off glia reservoir    

% Epileptic node
%--------------------------------------------------------------------------
epi     = 0;         	% Switches special epileptic node on or of
e       = fix(Ns / 2);  % Index of the abnormal node

eEE  = E.EE;            eEI  = 30;
eII  = E.II;            eIE  = E.IE;
eP   = E.P;           	eQ   = E.Q; 

% ODE options
%--------------------------------------------------------------------------
options     = odeset('MaxStep',0.5);
t_range     = 1:.005:100;
x_ini       = rand(Ns * 3,1) / 10;  % 0 -> near lower fix point; 


% Run simulations across different parameter values
%==========================================================================
clear B
oldx    = x_ini;

for d   = 1:dir
    
    if d == 1,      thisPR = Prange;        disp('Forward simulation'); 
    elseif  d == 2, thisPR = flip(Prange);  disp('Backward simulation'); 
    end
    
    i = 0; 
    
    for X = thisPR
        
        switch Pvary
            case 'EE',      EE  = X; 
            case 'EI',      EI  = X;
            case 'IE',      IE  = X; 
            case 'P',       P   = X;
            case 'c_EG',    EG  = X;
        end

        % Convert parameter values for use in model function
        %==================================================================
        % Excitatory connections within and between sources
        %------------------------------------------------------------------
        if connected    

            % Creates connectivity matrix A
            %--------------------------------------------------------------
            xloc        = 1:Ns;
            D           = dist(xloc);
            mu          = 0;
            sg          = 5;
            rdfact      = 0.0;      % Scales additional random connectivity
            A           = glia_dist2cnx(D, mu, sg, rdfact, Ns);
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
        
        c_EG    = zeros(Ns); 
        for c = 1:length(c_EG),     c_EG(c,c) = EG;     end
        T.c_EG  = c_EG;
        
        c_GE    = zeros(Ns); 
        for c = 1:length(c_GE),     c_GE(c,c) = GE;     end
        T.c_GE  = c_GE;

        % Other parameters
        %------------------------------------------------------------------
        T.A_S       = A_S;
        T.P         = ones(Ns,1) * P;
        T.Q         = ones(Ns,1) * Q;

        % Special 'Epileptic' Node
        %------------------------------------------------------------------
        if epi
            
        if connected,   T.c_EE(e,e) = eEE * [ Ns * EE / sum(sum(A)) ];
%                         T.c_EE(:,e) = 1.5 * T.c_EE(:,e);    % Rescale connectivity from epileptic node
%                         T.c_EE(e,:) = 0.5 * T.c_EE(e,:);    % Rescale connectivity to epileptic node
        else            T.c_EE(e,e) = eEE; end
        T.c_EI(e,e)     = eEI;      
        T.c_II(e,e)     = eII;
        T.c_IE(e,e)     = eIE;
        
        end
        
        % Glia as syncytium (rather than single pop per region)
        %------------------------------------------------------------------
        T.sync = sync;
        
        % Run actual simulations
        %------------------------------------------------------------------
        i = i + 1;
        disp(['Step ' num2str(i) ' of ' num2str(length(Prange))]);

        [t,x]       = ode45(@glia_cortsheet, t_range, oldx, options, T);
        oldx        = x(end,:)';
        
        % Store in output structure
        %------------------------------------------------------------------
        B(i,d).x      = x;
        B(i,d).t      = t;
        B(i,d).T      = T;
        
    end
end
C{cond} = B; 
end

%% Plotting routines
%==========================================================================
% Bifurcation plot
%--------------------------------------------------------------------------
xNo     = 1:Ns;
xNo     = xNo * 3 - 1 ;

clear Flo Fhi Blo Bhi
for f = 1:length(B)
    halfway = ceil(size(B(f).x) / 2);
    
    Flo(f, 1:length(xNo)) = real(log(min(B(f,1).x(halfway:end,xNo))));
    Fhi(f, 1:length(xNo)) = real(log(max(B(f,1).x(halfway:end,xNo))));
    
    if size(B,2) == 2 
        Blo(f, 1:length(xNo)) = real(log(min(B(f,2).x(halfway:end,xNo))));
        Bhi(f, 1:length(xNo)) = real(log(max(B(f,2).x(halfway:end,xNo)))); 
    end
    
end

clf 

cols = flip(cbrewer('div', 'Spectral', length(xNo)));
for i = 1:length(xNo)
%     if exist('Blo'), subplot(3,1,1); end
    subplot(3,1,1)
    scatter(Prange, Flo(:,i), 50, cols(i,:), 'filled'); hold on
    scatter(Prange, Fhi(:,i), 50, cols(i,:), 'filled');
    title(['Forward simulation, parameter ' Pvary]);
    xlim([min(Prange)-0.5 max(Prange)+0.5]);
    
    if exist('Blo') 
        subplot(3,1,2);
        thisPR  = flip(Prange);
        scatter(thisPR, Blo(:,i), 50, cols(i,:), 'filled'); hold on
        scatter(thisPR, Bhi(:,i), 50, cols(i,:), 'filled');
        title(['Backward simulation, parameter ' Pvary]);
        xlim([-.5 max(Prange)+0.5]);
    end
    
    subplot(3,1,3)
    clear cs
    for p = 1:length(Prange)
        cs(p) = mean(mean(cov(B(p).x(:,3*[1:Ns]-2))));
    end
    scatter(Prange, cs); 
end

%% Time series plot
%--------------------------------------------------------------------------
clf
ploti   = 100;
Nruns   = size(B,1);
Ei      = 3*[1:Ns] - 1;

subplot(1,7,[1:6])

for e = 1:length(Ei)
    ee = Ei(e);
    plot(B(ploti,1).x(:,ee) + e*0.5); 
  	hold on
end
title(['Parameter value: ' num2str(Prange(ploti))]);

subplot(1,7,7)
Bcor = corr(B(ploti,1).x(:,Ei));
imagesc(Bcor, [0 1]);
axis square

%% State space plot
%--------------------------------------------------------------------------
B = C{1};
for b = 1:length(B)
    xl = size(B(b).x,1); xr = fix(xl/2):xl; 
    plot(B(b).x(xr,1), B(b).x(xr,2)), hold on
end

%%
c = 1;
d = 1;
cols = flip(cbrewer('div', 'Spectral', length(C{c}))); 
for  b = 1:length(C{c})
xl = size(C{c}(b,1).x,1); xr = fix(xl/2):xl; 
if d == 1,  bvals = ones(length(xr),1)*b;
else,       bvals = ones(length(xr),1) * (length(C{c})-b+1);   
end

xe = C{c}(b,d).x(xr,10);
xi = C{c}(b,d).x(xr,11);
plot3(bvals, xe, xi, 'color', cols(b,:)); hold on
% ylim([0 0.25]);
end