function serial_dist_discr = discretise_serial_distrib(f_serial_cheb,t_discr)
    
    % Discretise the serial interval distribution, f_serial_cheb (which
    % should be a Chebfun), to calculate the probability mass function of
    % the (whole number) difference between the dates of onset of the
    % source and recipient, using the method in Cori et al., "﻿A new
    % framework and software to estimate time-varying reproduction numbers
    % during epidemics", Am. J. Epidemiol. (2013). The discretised
    % distribution is evaluated at the times in t_discr.
    
    F_serial = cumsum(f_serial_cheb);
    u = chebfun('u',f_serial_cheb.domain);
    g = cumsum(u*f_serial_cheb);
    
    p1 = (t_discr+1).*F_serial(t_discr+1);
    p2 = -2*t_discr.*F_serial(t_discr);
    p3 = (t_discr-1).*F_serial(t_discr-1);
    p4 = 2*g(t_discr)-g(t_discr-1)-g(t_discr+1);
    
    p_discr = p1+p2+p3+p4;
    
    serial_dist_discr.t = t_discr;
    serial_dist_discr.p = p_discr;
end