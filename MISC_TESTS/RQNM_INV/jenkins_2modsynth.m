clc
clear all
addpath('../ROUTINES/')
addpath('../ROUTINES/HARMONIC/')
addpath('../ROUTINES/SOLVERS/')

set(0,'defaultAxesTickLabelInterpreter', 'default');
set(0,'defaultTextInterpreter','latex'); 
set(0, 'DefaultLegendInterpreter', 'latex'); 
set(0,'defaultAxesFontSize',16)

cfg = 2;
load(sprintf('./DATA/Jenkins_RQNM_cfg%d.mat', cfg), 'As', 'Cs', 'Omegas', 'Phi', 'GM', 'Qs', 'UlC', 'Knl', 'Cnl', 'Q1s', 'Q2s');

Nq = length(Qs);
Q1s = squeeze(Q1s);
Q2s = squeeze(Q2s);

% Append zero case
Qs = [0; Qs];
Q1s = [zeros(1, Nq); Q1s];
Q1s = [Q1s(:,1) Q1s];

Q2s = [zeros(Nq, 1), Q2s];
Q2s = [Q2s(1, :); Q2s];

tmp = Knl;
Knl = zeros(3, Nq+1, Nq+1);
Knl(1, 2:end, 2:end) = tmp(1, :, :);
Knl(1, 1, 1:end) = Omegas(1, 1)^2;
Knl(1, 2:end, 1) = Omegas(:, 1).^2;

Knl(2, 2:end, 2:end) = tmp(2, :, :);
Knl(2, 1, 1:end) = 0;
Knl(2, 1:end, 1) = 0;

Knl(3, 2:end, 2:end) = tmp(3, :, :);
Knl(3, 1:end, 1) = Omegas(1, 2)^2;
Knl(3, 1, 2:end) = Omegas(:, 2).^2;

tmp = Cnl;
Cnl = zeros(3, Nq+1, Nq+1);
Cnl(1, 2:end, 2:end) = tmp(1, :, :);
Cnl(1, 1, 1:end) = Cs(1, 1);
Cnl(1, 2:end, 1) = Cs(:, 1);

Cnl(2, 2:end, 2:end) = tmp(2, :, :);
Cnl(2, 1, 1:end) = 0;
Cnl(2, 1:end, 1) = 0;

Cnl(3, 2:end, 2:end) = tmp(3, :, :);
Cnl(3, 1:end, 1) = Cs(1, 2);
Cnl(3, 1, 2:end) = Cs(:, 2);

Phi{1} = [Phi{1}(1,:); Phi{1}];
Phi{2} = [Phi{2}(1,:); Phi{2}];

m12 = Phi{1}*GM.M*Phi{2}';

F = GM.L(end,:)';
frc = [Phi{1}*F Phi{2}*F]';
frc = [frc; frc*0];
NLM = struct('Qs', Qs, 'Knl', Knl, 'Cnl', Cnl, 'm12', m12);

%% Single Mode Synthesis
FEX = GM.L(end,:)';

Omsq_1 = @(famp) sqrt((Omegas(:, 1).^2-Cs(:,1).^2/2) + sqrt(Cs(:,1).^4/4-(Omegas(:,1).*Cs(:,1)).^2 + ((Phi{1}(2:end,:)*FEX*famp)./Qs(2:end)).^2));
Omsq_2 = @(famp) sqrt((Omegas(:, 1).^2-Cs(:,1).^2/2) - sqrt(Cs(:,1).^4/4-(Omegas(:,1).*Cs(:,1)).^2 + ((Phi{1}(2:end,:)*FEX*famp)./Qs(2:end)).^2));
AmpSc = sqrt(0.5)*abs(GM.NLTs.L*Phi{1}(2:end,:)');
Phs = rad2deg(angle(Phi{1}(2:end,:)*FEX));

%% Check Gradients
rng(1);
uw0 = rand(5,1);
uw0(1:4) = uw0(1:4)*100;
[R0, dRdU, dRdw] = TWOMDRESFUN(uw0, frc, NLM);

Jnum = zeros(4,4);
hv = zeros(5, 1);
hm = 1e-6;
for hi=1:4
    hv(hi) = 1;
    Rp = TWOMDRESFUN(uw0+hv*hm, frc, NLM);
    Rm = TWOMDRESFUN(uw0-hv*hm, frc, NLM);
   
    Jnum(:, hi) = (Rp-Rm)/(2*hm);
    hv(hi) = 0;
end

hv(5) = 1;
Rp = TWOMDRESFUN(uw0+hv*hm, frc, NLM);
Rm = TWOMDRESFUN(uw0-hv*hm, frc, NLM);
Jwnum = (Rp-Rm)/(2*hm);
hv(5) = 0;

disp('Analytical')
disp(table([dRdU dRdw]))
disp('Numerical')
disp(table([Jnum Jwnum]))

% %% Solving
% u0 = [1;1;1;1];
% w0 = 7000;
% 
% opt = optimoptions('fsolve', 'SpecifyObjectiveGradient', true, 'Display', 'iter');
% fsolve(@(u) TWOMDRESFUN([u; w0], frc, NLM), u0, opt)

% %% solve & continue
% Sopt = struct('jac', 'none');
% USs = solve_and_continue(u0, @(uw) TWOMDRESFUN(uw, frc, NLM), 7500, 9000, 1, Sopt);

%% Fresp
load(sprintf('./DATA/Jenkins_FRESP_cfg%d.mat', cfg), 'UCs', 'Fl', 'Fas', 'Wst', 'Wen', 'h', 'Nhc');

Fas = [1e5 5e5 10e5 20e5 35e5];
Sopt = struct('jac', 'none');

u0 = [1;1;1;1];

figure(10)
clf()
figure(20)
clf()
for fi=1:5
    % 2 mode synthesis
    Us = solve_and_continue(u0, @(uw) TWOMDRESFUN(uw, Fas(fi)*frc, NLM), Wst, Wen, 0.5, Sopt);

    % 1 mode synthesis
    om1 = Omsq_1(Fas(fi));  
    om2 = Omsq_2(Fas(fi));
    
    ir1 = find(imag(om1)==0); 
    ir2 = find(imag(om2)==0);
    ir2 = ir2(end:-1:1);
    
    % Post process 2 mode synthesis
    Q1s = sqrt(Us(1, :).^2+Us(3,:).^2);
    th1s = atan2(Us(1, :), Us(3, :));
    Q2s = sqrt(Us(2, :).^2+Us(4,:).^2);
    th2s = atan2(Us(2, :), Us(4, :));

    Ws = Us(end,:);
    
    Res = [interp1(Qs, GM.NLTs.L*Phi{1}', Q1s).*Us(1,:) + interp1(Qs, GM.NLTs.L*Phi{2}', Q2s).*Us(2,:);
        interp1(Qs, GM.NLTs.L*Phi{1}', Q1s).*Us(3,:) + interp1(Qs, GM.NLTs.L*Phi{2}', Q2s).*Us(4,:)];

    figure(10)
    semilogy([om1(ir1); om2(ir2)], [Qs(ir1)' Qs(ir2)']/sqrt(2), 'k--', 'LineWidth', 1); hold on
    semilogy(Ws, Q1s/sqrt(2), 'g-', 'LineWidth', 2); hold on
    semilogy(Ws, Q2s/sqrt(2), 'm-', 'LineWidth', 2); hold on
    
    figure(20)
    plot(UCs{fi}(end,:), ...
        sqrt(sum((kron(diag([1 sqrt(0.5)*ones(1,Nhc-1)]), GM.NLTs.L)*UCs{fi}(1:end-1, :)).^2)), ...
        'b-', 'LineWidth', 2); hold on
    plot(Ws, sqrt([0.5 0.5]*Res.^2), 'r-', 'lineWidth', 1.5); hold on
    

    
    plot([om1(ir1); om2(ir2)], [AmpSc(ir1).*Qs(ir1)' AmpSc(ir2).*Qs(ir2)'], 'k--', 'LineWidth', 1);    
end
figure(10)
xlim([Wst Wen])

figure(20)
xlim([Wst Wen])

%% Ornament Figures
figure(10)
xlabel('Frequency (rad/s)')
ylabel('Harmonic Modal Amplitude')

text(8700, 4e-1, 'Mode 1', 'fontsize', 16)
text(8700, 8e-4, 'Mode 2', 'fontsize', 16)

annotation('arrow', ([8000 7900]-Wst)/(Wen-Wst), [0.5 0.8])
annotation('arrow', ([7950 7880]-Wst)/(Wen-Wst), [0.2 0.4])

text(7875, 1.75e-2, '$||f_{ex}||$', 'fontsize', 16)
text(7800, 2e-4, '$||f_{ex}||$', 'fontsize', 16)

print(sprintf('./FIGS/MODALRESP_SYNTH2_cfg%d.eps', cfg), '-depsc')

figure(20)
xlabel('Frequency (rad/s)')
ylabel('RMS Amplitude (m)')

annotation('arrow', ([8000 7900]-Wst)/(Wen-Wst), [0.45 0.9]/3)
text(7750, 0.7, '$||f_{ex}||$', 'fontsize', 16)
print(sprintf('./FIGS/FRESP_SYNTH2_cfg%d.eps', cfg), '-depsc')