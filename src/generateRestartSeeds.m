function samples = generateRestartSeeds(m,C,N)
    % Generate deterministic Gaussian restart seeds via a Sobol sequence.
    % Produces space-filling samples concentrated within a 3-sigma region.
    %
    % Inputs:
    %   m - Mean vector of the target Gaussian.
    %   C - Covariance matrix of the target Gaussian.
    %   N - Number of restart samples.
    %
    % Outputs:
    %   samples - Matrix of size [numel(m) x N] containing restart points.

    M = numel(m);
    Ssob = sobolset(M,'Skip',N);
    % Draw samples from the Sobol set, which are uniformly distributed on
    % [0,1]
    U = net(Ssob, N); 
    % Rescale U to concentrate samples within the 3-sigma region
    pmin = normcdf(-3);              
    pmax = normcdf( 3);         
    U = pmin + (pmax - pmin) * U; 

    % Map Sobol points U(0,1) -> N(0,1)
    Z = norminv(U.');    
    % Map N(0,1) -> N(bp,Cb)
    samples = m + chol(C,'lower') * Z;
end