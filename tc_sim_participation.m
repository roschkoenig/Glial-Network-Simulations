function Pc = tc_sim_participation(A, mod)
    Nclust = max(mod);
    for a = 1:size(A,1)
        thisa = squeeze(A(a,:,:));
        imagesc(thisa);
        thisa = thisa - diag(diag(thisa));
        
        for n = 1:size(thisa,1)
            modnds  = find(mod(mod == mod(n)));
            othnds  = find(mod(mod ~= mod(n))); 
            wi_cnx  = sum(sum(thisa(n,modnds))) / length(modnds); 
            bt_cnx  = sum(sum(thisa(n,othnds))) / length(othnds);
            Pc(a,n) = bt_cnx / wi_cnx;
        end
    end

end