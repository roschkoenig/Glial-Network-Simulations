function L = tc_layout(Nside)


L(1).name   = 'SEEG few';
Nshanks     = ceil(sqrt(Nside));
[L(1).chans, L(1).label] = tc_seeg_layout(Nside, Nshanks);

L(2).name   = 'SEEG many';
Nshanks     = ceil(Nside / 2);
[L(2).chans, L(2).label] = tc_seeg_layout(Nside, Nshanks);

L(3).name   = 'Grid';
Ngrid       = ceil(Nside / 3);
L(3).chans  = tc_grid_layout(Nside, Ngrid); 

end

function [chans label] = tc_seeg_layout(Nside, Nshanks)

    chans = zeros(Nside);
    label = cell(Nside);
    lobound     = floor(min([5, Nside/2]));
    hibound     = floor(max([lobound, Nside/2]));
    xpos        = randsample(Nside, Nshanks); 

    for n = 1:Nshanks
        Nchans      = fix(range([lobound hibound]) * rand(1)) + lobound; 
        leftovers   = Nside - Nchans;
        starty      = randsample(leftovers+1,1);
        stopy       = starty + Nchans-1;
        ypos        = starty:stopy;
        chans(xpos(n), ypos) = 1;
        for y = 1:length(ypos)
            label{xpos(n), ypos(y)} = [char(64+xpos(n)) num2str(ypos(y))];
        end
    end
end

function chans = tc_grid_layout(Nside, Ngrid)
    chans   = zeros(Nside);
    pos     = randsample(Nside+1 - Ngrid, 2, 1); 
    chans([1:Ngrid] + pos(1), [1:Ngrid] + pos(2)) = 1;    
end