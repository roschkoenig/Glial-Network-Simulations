clear F
fs          = filesep;
F.base      = '/Users/roschkoenig/Dropbox/Research/1909 CB Dispatch Sims';
F.scripts   = [F.base fs 'Network Sims'];

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

b = 93; % 92
xl = size(C{1}(b,1).x,1); xr = fix(xl-5000):2:xl;
cols = (cbrewer('seq', 'YlOrRd', length(xr)));

e = C{1}(b,1).x(xr,10);
i = C{1}(b,1).x(xr,11);
g = C{1}(b,1).x(xr,1);
colormap(cols)
scatter3(e,i,g, 20, g, 'filled')
view([-30, 22]); 

