function [] = mainrun(varargin)
    addpath('../ROUTINES/')
    addpath('../ROUTINES/FEM/')
    addpath('../ROUTINES/HARMONIC/')
    addpath('../ROUTINES/CONTACTMODELS/')
    addpath('../ROUTINES/SOLVERS/')
    addpath('../ROUTINES/QUADRATURE/')

    %% Parameters for run: 'kt', 'kn', 'mu', 'gap'
    parprops = [struct('par', 'kt', 'quadfun', @(n) GPHWT(n), 'map', @(x) 10^(x+12), 'mu', 12, 'sig', 1); 
        struct('par', 'kn', 'quadfun', @(n) GPHWT(n), 'map', @(x) 10^(x+12), 'mu', 12, 'sig', 1);
        struct('par', 'mu', 'quadfun', @(n) LAGWT(n), 'map', @(x) x*0.25, 'mu', 0.25, 'sig', nan);
        struct('par', 'gap', 'quadfun', @(n) LAGWT(n), 'map', @(x) x*1e-4, 'mu', 1e-4, 'sig', 1e-4)];
    
    Nq_pces = [1 1 1 1];
    Ixs = [1 1 1 1];
    nxi = 1;
    pref = '';
    
    %% Input Arguments
    if nargin~=0
       nxi     = varargin{1}; 
       Nq_pces = repmat(varargin{2}, [1 4]);
       
       Ixs    = str2num(dec2base(nxi-1, Nq_pces(1), 4)')'+1;
       
       pref = './DATA/RUNS/QUAD';
    end

    %% Run
    RQNM_LPF3_PCEFUN(Ixs, nxi, parprops, Nq_pces, pref);
end