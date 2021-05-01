function f_gen = get_gen_dist_ferretti(f_inc,b_cond)
    
    % Calculate the generation time distribution for the Ferretti model, by
    % averaging the conditional infectiousness, b_cond, over the incubation
    % period distribution, f_inc.
    
    % This function requires the Chebfun package to execute (freely
    % available at https://www.chebfun.org/download/).
    
    t_inc1_range = [0,25,50,500];
    t_gen_range = [0,25,50];
    
    f1 = @(t_gen)integral(@(t_inc1)b_cond(t_gen-t_inc1,t_inc1).*f_inc(t_inc1),t_inc1_range(1),t_inc1_range(2),'AbsTol',1e-12,'RelTol',1e-10);
    f2 = @(t_gen)integral(@(t_inc1)b_cond(t_gen-t_inc1,t_inc1).*f_inc(t_inc1),t_inc1_range(2),t_inc1_range(3),'AbsTol',1e-12,'RelTol',1e-10);
    f3 = @(t_gen)integral(@(t_inc1)b_cond(t_gen-t_inc1,t_inc1).*f_inc(t_inc1),t_inc1_range(3),t_inc1_range(4),'AbsTol',1e-12,'RelTol',1e-10);
    
    f_gen_fun = @(t_gen_vec)arrayfun(@(t_gen)f1(t_gen)+f2(t_gen)+f3(t_gen),t_gen_vec);
    f_gen1 = chebfun(f_gen_fun,t_gen_range);
    f_gen = f_gen1/sum(f_gen1);
end