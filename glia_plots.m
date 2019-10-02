clear F
fs          = filesep;
F.base      = '/Users/roschkoenig/Dropbox/Research/1909 CB Dispatch Sims';
F.scripts   = [F.base fs 'Network Sims'];
addpath(genpath(F.base)); 
load([F.base fs 'Sims.mat']); 

%% State space plot
%--------------------------------------------------------------------------
B = C{1};
for b = 1:length(B)
    xl = size(B(b).x,1); xr = fix(xl/2):xl; 
    plot(B(b).x(xr,1), B(b).x(xr,2)), hold on
end

%% Plot representations of bifurcations
%==========================================================================
do3d = 1;
cols = flip(cbrewer('div', 'Spectral', length(C{1}))); 
pval = linspace(-80,0,length(C{1}));

for c = 1:2, figure(c), 
for d = 1:2,
for b = 1:1:length(C{c})

xl = size(C{c}(b,1).x,1); xr = fix(xl-3000):2:xl-2000; 
if d == 1,  bvals = ones(length(xr),1)*pval(b);
else,       bvals = ones(length(xr),1) * (pval(end-b+1));   
end

xe = C{c}(b,d).x(xr,10);
xi = C{c}(b,d).x(xr,11);
    
% 2D Bifurcation plot
%--------------------------------------------------------------------------
if ~do3d
    subplot(2,1,d)
    scatter(bvals(1), log(max(xe)), 'filled', 'k'), hold on
    scatter(bvals(1), log(min(xe)), 'filled', 'k')
    ylim([-3.5 0.5]); 
end

% 3D Bifurcation plot
%--------------------------------------------------------------------------
if do3d && d == 1
    if c == 1, cols = (cbrewer('seq', 'Oranges', length(C{c}))); 
    else,      cols = (cbrewer('seq', 'Purples', length(C{c}))); 
    end
    plot3(bvals, xe, xi, 'color', cols(b,:)); hold on
    view([325 23]); 
    ylim([0 0.75]);
    zlim([0 1]); 
    xlim([-60 0])
    
    grid on
end
end
end
end

%% Plot example time series
%--------------------------------------------------------------------------
exps = [40, 90, 94]; 
for ex = 1:length(exps),    subplot(length(exps), 1, ex)
    c   = 1;
    b   = exps(ex); 
    xl  = size(C{1}(b,1).x,1); xr = fix(xl-15000):2:xl;

    eid = [1:size(C{1}(1).x,2)/3] * 3 - 2; 
    e   = C{c}(b,1).x(xr,eid(1));
    g   = C{c}(b,1).x(xr,eid(1)+2); 
    for k = 1:size(e,2), plot(e(:,k)-mean(e(:,k)) + k/3, 'color', 'k'), hold on; end
    for k = 1:size(e,2), plot(g(:,k)-mean(g(:,k)) + k/3, 'color', 'r'), hold on; end
end

%% Plot simulations with different temporal profiles
%--------------------------------------------------------------------------
exps = 95; %[50 80 94 118]
clear S
for ex = 1:length(exps)    
    ex
    % Set up initial conditions and parameters
    %----------------------------------------------------------------------
    c       = 1;
    B       = C{c}(exps(ex),1); 
    xini    = B.x(end,:); 
    T       = B.T; 
    T.t     = B.t; 
    T.A_S   = 0.1; 
    
    % Define time varying parameter values
    %----------------------------------------------------------------------
    cond = 'ramp';   % ping, ramp
    
    switch cond
        case 'ping'
            on = 9001;     off = 9500;
            P_t          = ones(1,length(T.t)) * -.1;
            P_t(on:off)  = ones(1,off-on+1) * 1;
            
        case 'ramp' 
            on = 3001;     off = 12000;
            rmp         = zeros(1,length(T.t)); 
            rmp(on:off) = (0:[off-on])/[off-on];
            P_t         = ones(1,length(T.t)) * -0.5; 
            P_t         = P_t + 1.5*rmp; 
            P_t(off:off+3000) = ones(1,length(off:off+3000));
            
            R_t         = 20*rmp;
    end  
    T.P_t        = P_t;
%     T.R_t        = R_t; 
%     T.R          = 0; 
    
    options     = odeset('MaxStep',0.5);
    [t,x]       = ode45(@glia_cortsheet, T.t, xini, options, T);
    S(ex).t     = t;
    S(ex).x     = x;
    S(ex).T     = T;
end

for ex = 1:length(S)
    x = S(ex).x; 
    t = S(ex).t;
    T = S(ex).T; 
    
    cols = cbrewer('qual', 'Set1', 3); 
    clvec       = zeros(length(x(:,1)),3); 
    for i = 1:on-1, clvec(i,:) = cols(2,:); end
    for i = on:off, clvec(i,:) = cols(1,:); end
    for i = off+1:length(x(:,1)), clvec(i,:) = cols(3,:); end
    
    plr         = 500:length(x(:,1)); 
    subplot(length(exps), 3, [ex*3-2 ex*3-1]); 
        plot(x(plr,10)); 
        
    subplot(length(exps), 3, ex*3)
        scatter3(x(plr,10), x(plr,11), x(plr,3), [], clvec(plr,:), 'filled')
end

%%
s = 1;
figure(1)
xr = 1000:length(S(s).x(:,10)); 
subplot(4,1,1), plot(S(s).x(xr,10)); 
subplot(4,1,2), plot(S(s).x(xr,11)); 
subplot(4,1,3), plot(S(s).x(xr,3)); 
subplot(4,1,4), plot(P_t(xr)); 

figure(2)

cols = (cbrewer('seq', 'YlOrRd', length(S(s).x(:,3))));
if c == 2, gc = ones(length(S(s).x(:,1))); 
else,      gc = floor(   (S(s).x(:,3)-min(S(s).x(:,3)))/ ....
                     (max(S(s).x(:,3))-min(S(s).x(:,3))) * ...
                     (length(S(s).x(:,3))-1) )+1;
end

xr = 1000:1:length(S(s).x(:,1)); 
for xi = 1:length(xr)-1
    plot(S(s).x(xr(xi:xi+1),10), S(s).x(xr(xi:xi+1),11), ...
         'color', cols(gc(xi),:)); hold on
end
ylim([0 0.3]); xlim([0.02 0.22]);
% view(44,40)
% axis tight 
% grid on

%% Plot detailed example
xr = 4000:3:length(S(2).x(:,1))-4000;
% plot3(S(2).x(xr,10), S(2).x(xr,11), -S(2).x(xr,3))
for s = 1:10
    plot(S(2).x(xr,(s-1)*3+2) + s/4), hold on
end
    plot(S(2).x(xr,3) - mean(S(2).x(xr,3)))



%%
cols = (cbrewer('seq', 'YlOrRd', length(xr)));

e = C{c}(b,1).x(xr,10);
i = C{c}(b,1).x(xr,11);
ef = C{1}(b,1).x(xr,1); 
g = C{c}(b,1).x(xr,3);
colormap(cols)
subplot(2,1,1)
    scatter3(e,i,g, 50, g, 'filled')
    view([-30, 22]); 
subplot(2,1,2)
    plot(zscore(g)); hold on
    plot(zscore(e))

