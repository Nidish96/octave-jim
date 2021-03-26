function [] = Ca_savedat(DOF, fsamp)
% clc
% clear all

Nein = 8
load(sprintf('./MATS/%dIN_MATS.mat', Nein), 'SensorLocs', 'RECOV', 'Vrbms', 'M');
Lrbms = null(full(Vrbms'*M));
clear Vrbms M

type = 'WGN';
% DOF = 'Z';

Famps = [1 100 500 1000];

Ts      = cell(3, 1);
Urecs   = cell(3, 1);
Udrecs  = cell(3, 1);
Uddrecs = cell(3, 1);
Exc     = cell(3, 1);

for i=1:length(Famps)
    if famp>1
        fname = sprintf('./DATA/%dIN_%sRESP_%s%d_samp%d.mat', Nein, type, DOF, famp, log2(fsamp));
    else
        fname = sprintf('./DATA/%dIN_%sRESP_%s%d_samp%d.mat', Nein, type, DOF, log10(famp), log2(fsamp));
    end
   load(fname, 'T', 'Urec', 'Udrec', 'Uddrec', 'fext', 'ldof')
%  load(sprintf('./DATA/%dIN_%sRESP_%s%d_samp%d.mat', Nein, type, DOF, Famps(i), log2(fsamp)), 'T', ...
%       'U', 'Ud', 'Udd', 'fext', 'ldof')
%   Urec = RECOV*Lrbms*U;
%   Udrec = RECOV*Lrbms*Ud;
%   Uddrec = RECOV*Lrbms*Udd;
%   clear U Ud Udd

  fsamp = 1.0/(T(2)-T(1));
  
  Ts{i}      = T;
  Urecs{i}   = Urec;
  Udrecs{i}  = Udrec;
  Uddrecs{i} = Uddrec;
  Exc{i}     = interp1(0:(1/fsamp):1, fext, Ts{i});

  fprintf('%d: %f %d\n', i, Famps(i), ldof);
  disp(max(abs(Urec(1:3,:)), [], 2))
end


save(sprintf('./DATA/%dIN_%sRESP_%sDOFEX_samp%d.mat', Nein, type, DOF, log2(fsamp)), 'Ts', ...
     'Urecs', 'Udrecs', 'Uddrecs', 'Exc', 'SensorLocs', 'fsamp', 'ldof', ...
     'Famps');
end
