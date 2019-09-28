% TC Network Analysis (Simulation)
%==========================================================================
% This script will run some network analyses on simulated neuronal activity
% recorded via virtual SEEG implantations

% Housekeeping 
%==========================================================================
D       = tc_housekeeping; 
fs      = filesep;
Fsim    = [D.Fscripts fs '1 - Network Simulations'];
impl    = 2;    % Type of implantation used (1: few SEEG, 2: many SEEG, 3: grid)

try load([Fsim fs 'EN25_randwalk.mat']);    end     % Simulated Data
try load([Fsim fs 'SEEGs.mat']); end                % Simulated SEEGs

%==========================================================================
% Produce simulated SEEG recordings *LONG TIME TO RUN* 
%==========================================================================
% This section simulates SEEG implantations (randomly) from the same
% cortical sheet and detects modules within the functional connectivitiy
% matrix used for the subsequent step - is only run if file isnt already
% available. Should be parallelised next time. 
%==========================================================================

rng(25)
if ~exist('S')
maxs = 1000;
for s = 1:maxs
    disp(['Running simulation ' num2str(s) ' of ' num2str(maxs)]);
    % Extract data from simulation
    %----------------------------------------------------------------------
    Nside = fix(sqrt(size(B.P,1)));
    Ei    = 1:2:size(B.x,2);

    % Generate SEEG layout
    %----------------------------------------------------------------------
    L     = tc_layout(Nside);
    allx  = B.x(:,Ei)'; 
    dat   = allx(find(L(impl).chans),:);
    dat(:,20000) = dat(:,19999);    % Fudge an extra sample to the end
    
    % Preprocess
    %----------------------------------------------------------------------
    dat     = ft_preproc_rereference(dat, 'all');
    dat     = ft_preproc_bandpassfilter(dat, 200, [1 80]);

    % Sliding window and community detection     
    %----------------------------------------------------------------------
    Fs = 200;                                                              % This should be done outside the loop
    A = tc_windowslide(dat, Fs, L(impl).label, 'pearson', 2*Fs , 2*Fs);

    k = 2;
    Acom    = tc_communities(A(2:end,:,:), k);
    Aprob   = tc_tpm(Acom);
    ProbComs = tc_communities(Aprob,k);
    labels  = L(impl).label(find(L(impl).chans));

    % Pack for saving
    %----------------------------------------------------------------------
    S(s).impl   = coms;
end

save([Fsim fs 'SEEGs.mat'], 'S'); 
end

%==========================================================================
% Evaluate modularity detection
%==========================================================================
% In this part of the code we are evaluating how well the modularity
% detection applied here identifies underlying structure within the
% implantation (i.e. spatial distance and connectedness within modules vs
% between modules) 
%==========================================================================

Nside   = fix(sqrt(size(B.P,1)));
A       = B.T.c_EE;
wi_dist = [];
bt_dist = [];
wi_cnx  = [];
bt_cnx  = [];

% Calculate within and between module connections and spatial distance
%--------------------------------------------------------------------------
for s = 1:length(S)
    
    % Extract positions of SEEG elecrodes
    %----------------------------------------------------------------------
    Nclust       = max(max(S(s).impl));
    [rpos,cpos]  = find(S(s).impl);
    pos          = find(S(s).impl);
    
    % Calculte euclidean distances
    %----------------------------------------------------------------------
    dists        = triu(dist([rpos'; cpos']));
    dists        = dists(find(dists));
    bt_dist      = [bt_dist, mean(dists)];
    
    % Extract connectivity weights
    %----------------------------------------------------------------------
	a            = A(pos,pos);
    a            = a - diag(diag(a));
    bt_cnx       = [bt_cnx, mean(mean(a))];
    
    for n = 1:Nclust
        
        % Extract positions of SEEG elecrodes (module)
        %------------------------------------------------------------------
        [rpos, cpos] = find(S(s).impl == n);
        pos          = find(S(s).impl == n);
        
        % Calculate euclideon distances (module)
        %------------------------------------------------------------------        
        dists        = triu(dist([rpos'; cpos']));
        dists        = dists(find(dists));
        if isempty(dists),  mdis = NaN; 
        else                mdis = mean(dists);     end
        wi_dist      = [wi_dist, mdis];
        
        % Extract connectivity weights (module)
        %------------------------------------------------------------------
        a               = A(pos,pos);
        a               = a - diag(diag(a));
        if isempty(a),  cnx = NaN; 
        else,           cnx = mean(a(find(a)));     end
        wi_cnx          = [wi_cnx, cnx];
        
    end
end

% Plot distributions of within vs between measures
%--------------------------------------------------------------------------
Dat(1).dat  = wi_dist;
Dat(1).norm = bt_dist;
Dat(1).name = 'Distance';
Dat(1).xlim = 'zerolim';

Dat(2).dat  = wi_cnx;
Dat(2).norm = bt_cnx;
Dat(2).name = 'Connectedness';

tc_raincloud(Dat);

%% ==========================================================================
% Subselection of nodes from module - criteria comparison * TAKES A WHILE *
%==========================================================================
% Here we are comparing different strategies to select 'representative'
% nodes within specific modules - the feature we are comparin against is
% spatial distance and connectedness to the 'epileptic' node
%==========================================================================

% Identify relevant traces from recording
%--------------------------------------------------------------------------
Nside = fix(sqrt(size(B.P,1)));
Ei    = 1:2:size(B.x,2);
dat   = B.x(:,Ei)';     dat(:,20000) = dat(:,19999);    % Fudge an extra sample to the end   
Fs    = 200;

% Sliding window across the whole recording
%--------------------------------------------------------------------------
clear labels;
for r = 1:Nside, for c = 1:Nside,   labels{r,c} = ['(' num2str(r) ',' num2str(c) ')']; end, end
A = tc_windowslide(dat, Fs, labels, 'pearson', 2*Fs , 2*Fs);

% Go through layout and identify subselection according to criteria
%--------------------------------------------------------------------------
clear N
textprogressbar('Finding Nodes: ');
for s = 1:length(S)
    
    % Identify clusters
    %----------------------------------------------------------------------
    Nclust      = max(max(S(s).impl));
    [rpos,cpos] = find(S(s).impl);
    pos         = find(S(s).impl);
    modi        = find(S(s).impl(:));
    mod         = S(s).impl(modi);
    a           = A(:,pos,pos);
    
    % Participation index of nodes
    %----------------------------------------------------------------------
    Pc          = tc_sim_participation(a, mod);
    mPc         = mean(Pc,1);
    
    % Eigenvector centrality
    %----------------------------------------------------------------------
    for n = 1:Nclust
        [rpos, cpos]    = find(S(s).impl == n);
        pos             = find(S(s).impl == n);
        a               = A(:,pos,pos);  
        m               = find(mod == n);
        
        [val Pcid]      = max(mPc(m));
        [val Dynid]     = max(range(Pc(:,m))); 
        eigname         = tc_eigmost(a, labels(pos));
        
        % Pack for exporting
        %------------------------------------------------------------------
        N(s).eig{n}    = find(strcmp(labels, eigname));
        N(s).par{n}    = pos(Pcid);
        N(s).dyn{n}    = pos(Dynid);

    end
    textprogressbar(s/length(S) * 100)
end
textprogressbar('Done');

% Calculate normative data for comparison (remains unused)
%--------------------------------------------------------------------------
mockarray = ones(Nside);
[rpos cpos] = find(mockarray);
dists   = dist([rpos'; cpos']);
dists   = max(dists(:)) - dists;
cnx     = B.T.c_EE - diag(diag(B.T.c_EE));

clear nd ed ec pd pc dd dc
for s = 1:length(N)
    eig = [N(s).eig{:}];    ed(s) = min(dists(1,eig));  ec(s) = max(cnx(1,eig));
    par = [N(s).par{:}];    pd(s) = min(dists(1,par));  pc(s) = max(cnx(1,par));
    dyn = [N(s).dyn{:}];    dd(s) = min(dists(1,dyn));  dc(s) = max(cnx(1,dyn));
    
    for r = 1:100
        channels = find(S(s).impl);
        sel         = randsample(channels, length(eig));
        nd(s,r)     = min(dists(1,sel));     nc(s,r)      = max(cnx(1,sel)); 
    end
end

%% Plot distributions of measures derived from different subsel criteria
%==========================================================================
clear X
X(1).dat    = {ec, pc, dc};
X(2).dat    = {ed, pd, dd};

cols    = flip(cbrewer('qual', 'Set1', 5));

for i = 1:length(X)
xdat    = X(i).dat;
for xx = 1:length(xdat)
    
    x = xdat{xx};
    % Plot data distribution
    %----------------------------------------------------------------------
    subplot(length(X), 1, i);
    [f,xi]  = ksdensity(x);
    plot(xi,f, 'color', cols(xx,:), 'linewidth', 1.5);     hold on
    md  = median(x);
    [val my] = min(abs(xi - md));
    my  = f(my);
    plot([md md], [0 my], ':', 'color', cols(xx,:), 'linewidth', 1.5);
end
end

legend({'Eigenvector', 'EV Median', 'Participation', 'P Median', 'Dynamic Participation', 'DP Median'})

%% ==========================================================================
% Plot example implantations as graphs 
%==========================================================================
% this takes a few virtual SEEG implantations and overlays them (and the
% modules detected) onto a graph of the connections between individual
% nodes
%==========================================================================
cutoff  = 95;
ct      = prctile(B.T.c_EE(:), cutoff);
hi_i    = find(B.T.c_EE > ct);
lo_i    = find(B.T.c_EE < ct);
wts     = zeros(size(B.T.c_EE));
wts(hi_i) = 1;

% Highlight implantation layout
%--------------------------------------------------------------------------
cmap = cbrewer('qual', 'Dark2', 20);
maxs = 4;
for s = 1:maxs
    subplot(1,maxs,s);
    gA      = graph(wts, 'omitselfloops');
    pA      = plot(gA, 'layout', 'force');
    axis square
    highlight(pA,1)
    
    for m = 1:max(S(s).impl(:))
        [rpos, cpos]    = find(S(s).impl == m);
        pos             = find(S(s).impl == m);
        a               = A(:,pos,pos);  
        
        for nn = 1:length(pos)
        for n = nn:length(pos)
            pedge = [pos(nn), pos(n)];
        end
        end
        
        highlight(pA,pos, 'NodeColor', cmap(m,:), 'Markersize', 3);
        
        eigname         = tc_eigmost(a, labels(pos));
        eig(m)          = find(strcmp(labels, eigname));
    end
 
end

