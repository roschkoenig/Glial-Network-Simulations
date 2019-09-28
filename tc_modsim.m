function A = tc_modsim 

Nshanks = 10;    % Number of Shanks
canv    = 25;   % Canvas size

rng(24);

% Generate electrode layout
%--------------------------------------------------------------------------
subplot(2,3,[1 4])
[S E] = tc_modsim_layout(Nshanks, canv);

% Calculate spatial distance matrix
%--------------------------------------------------------------------------
allpos      = vertcat(S.cpos);
epipos      = vertcat(E.cpos);
D           = dist([allpos; epipos]');

clear label
for a = 1:size(allpos,1)
    label{a} = ['X' num2str(fix(allpos(a,1))) 'Y' num2str(fix(allpos(a,2)))];
end

% Convolve spatial distance with Gaussian weights for connectivity
%--------------------------------------------------------------------------
mu      = 0;
sg      = 5;
rdfact  = 0.1;      % Scales additional random connectivity
A       = tc_dist2cnx(D, mu, sg, rdfact, length(allpos));
Avis    = A(1:length(allpos), 1:length(allpos));

% Run community detection (just for visualisation sake) 
%--------------------------------------------------------------------------
B = tc_communities(Avis,fix(Nshanks/2));

% Plot results
%--------------------------------------------------------------------------
subplot(2,3,[2 3 5 6])
    tc_sortplot(Avis, B, label);
    
subplot(2,3,[1 4])
    tc_modsim_commplot(S,B);


% Equip each node with a Wilson-Cowan type oscillator sitting at noisy
% resting state

end

function [S E] = tc_modsim_layout(Nshanks, canv)

    % Define spatial location of individual SEEG contacts
    %----------------------------------------------------------------------

    clear S E

    llim = [10 20]; 
    llist = floor(( llim(2) - llim(1) ) .* rand(Nshanks,1) + llim(1));

    cols = cbrewer('qual', 'Dark2', Nshanks);

    for n = 1:Nshanks

        % Randomise x-position (burr-hole position) 
        %------------------------------------------------------------------
        S(n).x  = 0;
        xpos    = rand(1) * canv;
        while min(abs([S.x] - xpos)) < 1
            xpos = rand(1) * (canv-1);
        end
        S(n).x  = xpos;

        % Randomise y-position (implantation depth)
        %------------------------------------------------------------------
        leftspace = canv - llist(n) - 2;
        ypos      = rand * leftspace + 1;
        S(n).y    = [ypos ypos + llist(n)];

        % Plotting routines
        %------------------------------------------------------------------
        S(n).col  = cols(n,:);

        % Plot shanks
        plot([S(n).x S(n).x], [S(n).y(1), S(n).y(2)], 'color', S(n).col, 'linewidth', 2); 
        hold on

        % Plot contacts
        xs = ones(llist(n) + 1,1) * S(n).x;
        ys = [0:llist(n)]' + S(n).y(1);

        S(n).cpos = [xs, ys];
        scatter(S(n).cpos(:,1), S(n).cpos(:,2), 25, S(n).col);

        xlim([0 canv]);
        ylim([0 canv]);

    end
    
    % Add hidden epileptogenic nodes
    %======================================================================
    oldpos = [0 0];
    for e = 1:2
        E(e).cpos = oldpos;
        while max(dist([E(e).cpos; oldpos]')) < fix(canv/2.5)
            E(e).cpos = rand(1,2) * canv;
        end
        oldpos = E(e).cpos;
        scatter(E(e).cpos(1), E(e).cpos(2), 70, 'k', 'filled');
    end

end

% function A = tc_modsim_dist2cnx(D, mu, sg, rf, Nc)
%     
%     mdist   = mean(mean(D));
%     A      	= zeros(size(D,1), size(D,2));      % Initialise matrix
%     A(:)    = normpdf(D(:), 0, mdist/3);        % Fill matrix with weights
%     
%     Enorm       = normpdf(D(Nc+1:end, :), 0, mdist);
%     Enorm       = max(max(A)) * Enorm / (max(max(Enorm)));
%     A(Nc+1:end, :) = Enorm;
%     
%     Enorm       = normpdf(D(:, Nc+1:end), 0, mdist);
%     Enorm       = max(max(A)) * Enorm / (max(max(Enorm)));
%     A(:, Nc+1:end) = Enorm;
%    
%     % Make symmetric noise matrix
%     %----------------------------------------------------------------------
%     R       = rf * rand(length(A));
%     for rr = 1:length(R)
%     for rc = 1:length(R)
%         R(rr,rc) = R(rc,rr);
%         if rr == rc
%             R(rr,rc) = 0;
%         end
%     end
%     end
%     
%     % Combine noise and adjacency matrices
%     %----------------------------------------------------------------------
%     A       = A + R;
%     
%     A       = A / max(max(A));                  % Normalise matrix
%     for a = 1:length(A)                         % Ensure that diagonal is = 1
%         A(a,a) = 1;
%     end
% end

function tc_modsim_commplot(S,B)

% This function plots the channel locations community-representative
% colours

i       = 0;
cols    = jet(max(B));

for s = 1:length(S)
for c = 1:length(S(s).cpos)
    i = i+1;
    hold on
    scatter(S(s).cpos(c,1), S(s).cpos(c,2), 50, cols(B(i),:), 'filled');
end
end


end