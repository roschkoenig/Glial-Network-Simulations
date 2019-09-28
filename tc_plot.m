function tc_plot(P)

    x   = P.x;
    if isfield(P, 'cmap'),      cmap    = P.cmap;   else,   cmap    = copper(200);      end
    if isfield(P, 'dist'),      d       = P.dist;   else,   d       = 1:size(P.x,2);    end
    if isfield(P, 'trange'),    trange  = P.trange; else,   trange  = 1:size(P.x,1);    end
    if isfield(P, 'yscale'),    sc      = P.yscale; else,   sc      = 1;                end
    if isfield(P, 'Npoints'),   Npoints = P.Npoints; else,  Npoints = length(trange);   end
    if isfield(P, 'plsel'),     plsel   = P.plsel;  else,   plsel   = 1;                end
    if isfield(P, 't'),         t       = P.t;      else    t = trange;                 end
        

    lwidth              = 0.6;
    scalefactor         = mean(range(x(trange,:))) / (size(x,2) * sc / plsel);
    if Npoints < length(trange),    tsubset = fix(linspace(trange(1), trange(end), Npoints));
    else,                           tsubset = trange;
    end

%     plot(x(tsubset,1) - mean(x(tsubset,6)), 'color', 'r', 'linewidth', lwidth); hold on
    for r = 1:size(x,2)
        if mod(r,plsel) == 0
            plx    = x(tsubset,r) - mean(x(tsubset,r));    % extract data and BL correct
            plot(t(tsubset), plx - r*scalefactor, 'color', cmap(d(r),:), 'linewidth', lwidth); hold on
        end
    end
    xlim([-Inf Inf]);
    ylim([-Inf Inf]);
%     axis off