clc
clear all
addpath('../ROUTINES')
addpath('../ROUTINES/FEM')
addpath('../ROUTINES/CONTACTMODELS')
addpath('../ROUTINES/QUASISTATIC')
addpath('../ROUTINES/TRANSIENT')
addpath('../ROUTINES/HARMONIC')
addpath('../ROUTINES/SOLVERS')

model = 'BRB';

Runs = {'0001', '001', '01', '2.5', '05', '10', '15', '25', '50'};
load('./DATS/TRANSPROC.mat', 'Amps', 'Freqs', 'Damps')

type = 1;  % FRF
% type = 0;  % Resp

% ir = 1;

aa = gobjects(size(Runs));
for ir=1:length(Runs)
    load(sprintf('./DATS/HBCONT_R%s.mat',Runs{ir}), 'UwC', 'dUdwC', 'R', 'fa', 'h')  % Load Data

    Nhc = uint32(sum(h==0)+2*sum(h~=0));

    Ws = UwC(end,:);
    Resp = kron(eye(Nhc), R(3,:))*UwC(1:end-1,:);
    dResp = kron(eye(Nhc), R(3,:))*dUdwC(1:end-1,:);

    Np = length(Ws);
    Nq = 5;
    Q = EBBM2D_ND2QP(diff(Ws), Nq);

    Ws_q = zeros(Np+(Np-1)*Nq,1);
    Resp_q = zeros(Nhc, Np+(Np-1)*Nq,1);
    vpts = 1:(Nq+1):Np+(Np-1)*Nq;
    ipts = setdiff(1:Np+(Np-1)*Nq, vpts);

    Ws_q(vpts) = Ws;
    Ws_q(ipts) = Q(1:2:end,1:3:end)*Ws';
    Resp_q(:, vpts) = Resp;
    Resp_q(:, ipts) = (Q(2:2:end,2:3:end)*Resp')' + (Q(2:2:end,3:3: ...
                                                      end)*dResp')';

    figure(1)
    if ir==1
        clf()
    end
    if type==1
        semilogy(Ws/2/pi, sqrt(sum(Resp(2:3,:).^2,1))/fa, 'k.'); hold on
        aa(ir) = plot(Ws_q/2/pi, sqrt(sum(Resp_q(2:3,:).^2,1))/fa, '-', ...
                      'LineWidth', 2); hold on    
        legend(aa(ir), sprintf('F = %.3f N', fa))
    else 
        semilogy(Ws/2/pi, sqrt(sum(Resp(2:3,:).^2,1)), 'k.'); hold on
        aa(ir) = plot(Ws_q/2/pi, sqrt(sum(Resp_q(2:3,:).^2,1)), '-', ...
                      'LineWidth', 2); hold on    
        legend(aa(ir), sprintf('F = %.3f N', fa))        
    end

    xlim([145 170])
    xlabel('Forcing Frequency (Hz)')
    ylabel('Response FRF (m/N)')

    figure(2)
    if ir==1
        clf()
    end
    plot(Ws/2/pi, rad2deg(angle([1 -1j]*Resp(2:3,:))), 'k.'); hold on
    plot(Ws_q/2/pi, rad2deg(angle([1 -1j]*Resp_q(2:3,:))), '-', 'LineWidth', 2); hold on

    xlim([145 170])
    xlabel('Forcing Frequency (Hz)')
    ylabel('Respone Phase (degs)')
end

if type~=1
    figure(1)
    
end

if type==1
    figure(1)
    legend(aa, 'Location', 'best')
    print('./FIGS/FRF_AMP.eps', '-depsc')
    
    figure(2)
    print('./FIGS/FRF_PHASE.eps', '-depsc')
else
    figure(1)
    bb = plot(Freqs, Amps./(2*pi*Freqs).^2, 'k--', 'LineWidth', 2);
    legend(bb, 'PFF');
    
    legend([aa bb], 'Location', 'best')
    print('./FIGS/FRESP_AMP.eps', '-depsc')
end