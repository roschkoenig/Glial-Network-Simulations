function tc_raincloud(D)
clf
gcf

try     cols = cbrewer('qual', 'Set1', length(D)); 
catch   cols = jet(length(D));  end
    
for d = 1:length(D)
    subplot(length(D), 1, d);
    x = D(d).dat;
    lilx(d)    = min(x);
    bigx(d)    = max(x);
    
    % Plot data distribution
    %----------------------------------------------------------------------
    [f,xi] = ksdensity(x);
    plot(xi,f, 'color', cols(d,:), 'linewidth', 1.5);     hold on

    % Plot data points
    %----------------------------------------------------------------------
    jitscale = range(f) / 3;
    jit = jitscale * rand(1,length(x));
    scatter(x, jit - jitscale*1.5, 'markerfacecolor', cols(d,:), 'markerfacealpha', .1, 'markeredgecolor', 'none');
    
    % Plot boxplot
    %----------------------------------------------------------------------
    p25     = prctile(x,25); 
    p75     = prctile(x,75); 
    p50     = prctile(x,50);
    p95     = prctile(x,95);
    p05     = prctile(x,5);
    
    loline  = -1.25*jitscale;
    hiline  = -.75*jitscale;
    tpline  = -.5*jitscale;
    btline  = -1.5*jitscale;
    mdline  = -1*jitscale;
    
    plot([p25, p25],[loline hiline], 'color', 'k')
    plot([p75, p75],[loline hiline], 'color', 'k')
    plot([p25, p75],[loline loline], 'color', 'k')
    plot([p25, p75],[hiline hiline], 'color', 'k')
    plot([p50, p50],[loline hiline], 'color', 'k', 'linewidth', 2);
    plot([p05, p95],[mdline mdline], 'color', 'k', 'linewidth', 1.5);
    
    % Plot normative distribution markers
    %----------------------------------------------------------------------
    if isfield(D, 'norm')
        mnorm = mean(D(d).norm);
        plot([mnorm mnorm], [btline tpline], ':', 'linewidth', 1.2, 'color', [.2 .2 .2]);
        [h p]   = ttest2(D(d).norm, D(d).dat);
        
        if h
            if p <.01,  siglab = '**';
            else,       siglab = '*'; end
            text(p50, tpline + .2*jitscale, siglab, 'HorizontalAlignment', 'center', 'FontSize', 20);
        end
            
    end
    
    % Plot settings
    %----------------------------------------------------------------------
    ylim([-.6*max(f) 1.1*max(f)])
    title(D(d).name);
end

if isfield(D, 'xlim')
    switch D(1).xlim
        case 'same',    lims = [min(lilx) max(bix)];
        case 'zerolim', lims = [0 Inf];
    end
else
    lims = [-Inf Inf];
end

for d = 1:length(D)
    subplot(length(D), 1, d)
    xlim(lims)
end

end