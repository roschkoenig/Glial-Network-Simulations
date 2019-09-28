function A = glia_dist2cnx(D, mu, sg, rf, Nc)
    
    mdist   = mean(mean(D));
    A      	= zeros(size(D,1), size(D,2));      % Initialise matrix
    A(:)    = normpdf(D(:), 0, mdist/3);        % Fill matrix with weights
    
    if Nc < length(A)
        Enorm       = normpdf(D(Nc+1:end, :), 0, mdist);
        Enorm       = max(max(A)) * Enorm / (max(max(Enorm)));
        A(Nc+1:end, :) = Enorm;

        Enorm       = normpdf(D(:, Nc+1:end), 0, mdist);
        Enorm       = max(max(A)) * Enorm / (max(max(Enorm)));
        A(:, Nc+1:end) = Enorm;
    end
   
    % Make symmetric noise matrix
    %----------------------------------------------------------------------
    R       = rf * rand(length(A));
    for rr = 1:length(R)
    for rc = 1:length(R)
        R(rr,rc) = R(rc,rr);
        if rr == rc
            R(rr,rc) = 0;
        end
    end
    end
    
    % Combine noise and adjacency matrices
    %----------------------------------------------------------------------
    A       = A + R;
    
    A       = A / max(max(A));                  % Normalise matrix
    for a = 1:length(A)                         % Ensure that diagonal is = 1
        A(a,a) = 1;
    end
end