D       = tc_housekeeping; 
fs      = filesep;
Fsim    = [D.Fscripts fs '1 - Network Simulations'];
load([Fsim fs 'EN25_top_right.mat']);
%%
x = ft_preproc_bandpassfilter(B.x, 1000, [1 120]);

% Set up colour legend
%--------------------------------------------------------------------------
Ei      = 1:2:size(x,2);
Sarray  = sqrt(length(Ei));
i       = 0;

for nr = 1:Sarray
for nc = 1:Sarray
    i   = i + 1;
    xloc(i,:) = [nr,nc];
end
end
D       = dist(xloc');

edi    = 1; %fix(Sarray / 2);
idx     = xloc == edi;
idx     = intersect(find(idx(:,1)), find(idx(:,2)));
dis     = D(idx,:);
dis     = (dis - min(dis)) / range(dis);    % Normalise 0 to 1
dis     = fix(99 * dis) + 1;                % Range 1 - 100;

%% Plot the connectivity matrix
%--------------------------------------------------------------------------
cmap = flip(cbrewer('seq', 'YlGnBu', 100));
colormap(cmap);
imagesc(log(B.T.c_EE), [-3 -2]);
axis square


%%
[sorted sorting]    = sort(dis);
Ex                  = x(:,Ei);
Ex                  = Ex(:,sorting);

selection           = 5;       % plot only ever nth electrode
trange              = 10:length(B.t); %5500:6500;
scalemore           = .1;
plotpoints          = 5000;

% Create plot object
%--------------------------------------------------------------------------
ts      = [10 length(B.t);  6000 8000];

P           = [];           P.t         = B.t;
P.x         = Ex;           P.plsel 	= 10;         
P.dist      = sorted;       P.cmap      = (bone(200));          
P.yscale    = .1;           P.Npoints   = 5000;

for t = 1:size(ts,1)
    P.trange = ts(t,1) : ts(t,2);
    subplot(size(ts,1),1,t)
        tc_plot(P);
end

%% 
bpE = x(:,Ei) - mean(x(:,Ei));
ee  = [10010 10011 10012 10013 10014];

for i = 1:length(ee)
    e       = ee(i);
    thisx   = reshape(bpE(e,:), Sarray, Sarray);
    subplot(1,length(ee), i)
    imagesc(thisx, [-.02 .02])
    axis square
    cmap = cbrewer('div', 'Spectral', 100);
    colormap(cmap)
end

%% 
cutoff  = 95;
ct      = prctile(B.T.c_EE(:), cutoff);
hi_i    = find(B.T.c_EE > ct);
lo_i    = find(B.T.c_EE < ct);
wts     = zeros(size(B.T.c_EE));
wts(hi_i) = 1;

gA      = graph(wts, 'omitselfloops');
pA      = plot(gA, 'layout', 'force3');
highlight(pA,1)




