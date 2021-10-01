    Cnl = zeros(2, 2, Nq, Nq);
    for q1i=1 % :Nq
        q1t = cos(Omegas(q1i, 1)*T)*Qs(q1i);
        q1dot = -Omegas(q1i, 1)*sin(Omegas(q1i, 1)*T)*Qs(q1i);
        for q2i=1:10 % Nq
            q2t = cos(Omegas(q2i, 2)*T)*Qs(q2i);
            q2dot = -Omegas(q2i, 2)*sin(Omegas(q2i, 2)*T)*Qs(q2i);
            
            ut = Phi{1}(q1i,:).*q1t + Phi{2}(q2i,:).*q2t;
            udot = Phi{1}(q1i, :).*q1dot + Phi{2}(q2i, :).*q2dot;
            
            % Initialize slider states
            GM.NLTs.up = 0;
            GM.NLTs.fp = 0;
            for ti=1:length(T)
                [force_t(ti,:), ~, ~, GM] = GM.NLFORCE(T(ti), ut(ti,:)', udot(ti,:)', T(ti)-T(2));  % Nonlinear force
                force_t(ti, :) = force_t(ti, :) + ...
                    udot(ti, :)*GM.C' + ut(ti, :)*GM.K'; % Adding Linear forcing dissipation
            end
            
            % Regression to fit c11, c12, c22 (symmetric modal C)
            fmodal_t = ([Phi{1}(q1i, :); Phi{2}(q2i, :)]*force_t')';  % (Nt, 2)
            fmodal_t = fmodal_t - [q1t q2t]*diag([Omegas(q1i, 1)^2; Omegas(q2i, 2)^2]);  % Removing "conservative part"

            reg = [q1dot q2dot zeros(length(T),1);zeros(length(T),1) q1dot q2dot]\[fmodal_t(:,1); fmodal_t(:,2)];
            Cnl(:, :, q1i, q2i) = [reg(1) reg(2);reg(2) reg(3)];
            
            fprintf('(%d, %d)\n', q1i, q2i);
        end
    end