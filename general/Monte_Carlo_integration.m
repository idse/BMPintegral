function I = Monte_Carlo_integration(f, m, L, Npoints, D)
    % Monte Carlo integration 
    % I = Monte_Carlo_integration(f, L, Npoints)
    % 
    % f : function handle 
    % m : center of domain
    % L : linear size of domain
    % D : dimension 
    % Npoints : number of points to sample the function at
    
    % transpose mean if necessary
    if size(m,1) ~= 1
        m = m';
    end
    
    x = m + (rand([Npoints D]) - 0.5)*2*L; 
    tic
    I = zeros([Npoints 1]);
    for i = 1:Npoints
        I(i) = f(x(i,:)');
    end
    I = nansum(I(:)); % unfortunately matlab gives 0 log 0 = NaN instead of 0 so use nansum
    toc
    
    I = (2*L)^D*I/Npoints;
end