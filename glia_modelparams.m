function E = glia_modelparams(name)
    switch name
        case 'hopf' 
            E.EE  = 24;           E.EI  = 40;
            E.II  = 0;            E.IE  = -20;
            E.P   = 1.5;          E.Q   = -2; 
            
        case 'snic'
            E.EE  = 23;           E.EI  = 35;
            E.II  = 0;            E.IE  = -15;
            E.P   = 0.5;          E.Q   = -5;             
    end  
end
