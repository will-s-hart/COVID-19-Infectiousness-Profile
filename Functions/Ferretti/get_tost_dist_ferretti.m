function f_tost = get_tost_dist_ferretti(f_inc,b_cond)
    
    % Calculate the TOST distribution for the Ferretti model, by averaging
    % the conditional infectiousness, b_cond, over the incubation period
    % distribution, f_inc.
    
    % This function requires the Chebfun package to execute (freely
    % available at https://www.chebfun.org/download/).

    t_inc1_range = [0,25,50,500];
    t_tost_range = [-50,-25,0,25,50];
    
    f1 = @(t_tost)integral(@(t_inc1)b_cond(t_tost,t_inc1).*f_inc(t_inc1),t_inc1_range(1),t_inc1_range(2),'AbsTol',1e-12,'RelTol',1e-10);
    f2 = @(t_tost)integral(@(t_inc1)b_cond(t_tost,t_inc1).*f_inc(t_inc1),t_inc1_range(2),t_inc1_range(3),'AbsTol',1e-12,'RelTol',1e-10);
    f3 = @(t_tost)integral(@(t_inc1)b_cond(t_tost,t_inc1).*f_inc(t_inc1),t_inc1_range(3),t_inc1_range(4),'AbsTol',1e-12,'RelTol',1e-10);
    
    f_tost_fun = @(t_tost_vec)arrayfun(@(t_tost)f1(t_tost)+f2(t_tost)+f3(t_tost),t_tost_vec);
    f_tost_inc1 = chebfun(f_tost_fun,t_tost_range);
    f_tost = f_tost_inc1/sum(f_tost_inc1);
end