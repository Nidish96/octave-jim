function [] = Ca_savedat(DOF)
% clc
% clear all

Nein = 8
load(sprintf('./MATS/%dIN_MATS.mat', Nein), 'SensorLocs', 'RECOV', 'Vrbms', 'M');
Lrbms = null(full(Vrbms'*M));
clear Vrbms M

type = 'WGN';
% DOF = 'Z';

Famps = [100 500 1000];

Ts      = cell(3, 1);
Urecs   = cell(3, 1);
Udrecs  = cell(3, 1);
Uddrecs = cell(3, 1);
Exc     = cell(3, 1);

for i=1:length(Famps)
%   load(sprintf('./DATA/%dIN_%sRESP_%s%d.mat', Nein, type, DOF, Famps(i)), 'T', ...
%        'Urec', 'Udrec', 'Uddrec', 'fext', 'ldof')
  load(sprintf('./DATA/%dIN_%sRESP_%s%d.mat', Nein, type, DOF, Famps(i)), 'T', ...
       'U', 'Ud', 'Udd', 'fext', 'ldof')
   Urec = RECOV*Lrbms*U;
   Udrec = RECOV*Lrbms*Ud;
   Uddrec = RECOV*Lrbms*Udd;
   clear U Ud Udd

  fsamp = 1.0/(T(2)-T(1));
  
  Ts{i}      = T;
  Urecs{i}   = Urec;
  Udrecs{i}  = Udrec;
  Uddrecs{i} = Uddrec;
  Exc{i}     = interp1(0:(1/fsamp):1, fext, Ts{i});

  fprintf('%d: %f %d\n', i, Famps(i), ldof);
  disp(max(abs(Urec(1:3,:)), [], 2))
end


save(sprintf('./DATA/%dIN_%sRESP_%sDOFEX.mat', Nein, type, DOF), 'Ts', ...
     'Urecs', 'Udrecs', 'Uddrecs', 'Exc', 'SensorLocs', 'fsamp', 'ldof', ...
     'Famps');
end
