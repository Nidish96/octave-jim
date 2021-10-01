function [R, dRdUw, dRdw] = TWOMDRESFUN(uw, frc, NLM)
    U = uw(1:4);  % [q1c; q2c; q1s; q2c]
    w = uw(5); 
    
    Q1 = sqrt(U(1)^2+U(3)^2);
    if Q1~=0
        dQ1 = [U(1) 0 U(3) 0]/Q1;
    else
        dQ1 = [0 0 0 0];
    end
    Q2 = sqrt(U(2)^2+U(4)^2);
    if Q2~=0
        dQ2 = [0 U(2) 0 U(4)]/Q2;
    else
        dQ2 = [0 0 0 0];
    end
    
    knl = zeros(3,1);
    cnl = zeros(3,1);
    m12 = 0;
    
    fm = zeros(4, 1);
    
    dknl = zeros(3,2);  % [ddQ1, ddQ2]
    dcnl = zeros(3,2);  % [ddQ1, ddQ2]
    dm12 = zeros(1,2);
    
    dfm = zeros(4, 2);
    
    % Interpolate
    Qs = NLM.Qs;
    q1i = find((Qs(1:end-1)-Q1).*(Qs(2:end)-Q1)<=0); q1i = q1i(1);
    if isempty(q1i)
        error('Q1 out of bounds')
    end
    q2i = find((Qs(1:end-1)-Q2).*(Qs(2:end)-Q2)<=0); q2i = q2i(1);
    if isempty(q2i)
        error('Q2 out of bounds')
    end
    
    N1 = [(Qs(q1i+1)-Q1) (Q1-Qs(q1i))]/(Qs(q1i+1)-Qs(q1i));  % Interpolants
    dN1= [-1 1]/(Qs(q1i+1)-Qs(q1i));
    N2 = [(Qs(q2i+1)-Q2) (Q2-Qs(q2i))]/(Qs(q2i+1)-Qs(q2i));  % Interpolants
    dN2= [-1 1]/(Qs(q2i+1)-Qs(q2i));
    
    N12 = [N1(1)*N2(1) N1(1)*N2(2) N1(2)*N2(1) N1(2)*N2(2)];  % [x1y1 x1y2 x2y1 x2y2]
    dNQ1= [dN1(1)*N2(1) dN1(1)*N2(2) dN1(2)*N2(1) dN1(2)*N2(2)];
    dNQ2= [N1(1)*dN2(1) N1(1)*dN2(2) N1(2)*dN2(1) N1(2)*dN2(2)];
    
    knl = NLM.Knl(:, q1i, q2i)*N12(1) + NLM.Knl(:, q1i, q2i+1)*N12(2) + ...
        NLM.Knl(:, q1i+1, q2i)*N12(3) + NLM.Knl(:, q1i+1, q2i+1)*N12(4);
    cnl = NLM.Cnl(:, q1i, q2i)*N12(1) + NLM.Cnl(:, q1i, q2i+1)*N12(2) + ...
        NLM.Cnl(:, q1i+1, q2i)*N12(3) + NLM.Cnl(:, q1i+1, q2i+1)*N12(4);
    m12 = NLM.m12(q1i, q2i)*N12(1) + NLM.m12(q1i, q2i+1)*N12(2) + ...
        NLM.m12(q1i+1, q2i)*N12(3) + NLM.m12(q1i+1, q2i+1)*N12(4);
    
    dknl = NLM.Knl(:, q1i, q2i).*[dNQ1(1) dNQ2(1)] + ...
        NLM.Knl(:, q1i, q2i+1)*[dNQ1(2) dNQ2(2)] + ...
        NLM.Knl(:, q1i+1, q2i).*[dNQ1(3) dNQ2(3)] + ...
        NLM.Knl(:, q1i+1, q2i+1)*[dNQ1(4) dNQ2(4)];
	dcnl = NLM.Cnl(:, q1i, q2i).*[dNQ1(1) dNQ2(1)] + ...
        NLM.Cnl(:, q1i, q2i+1)*[dNQ1(2) dNQ2(2)] + ...
        NLM.Cnl(:, q1i+1, q2i).*[dNQ1(3) dNQ2(3)] + ...
        NLM.Cnl(:, q1i+1, q2i+1)*[dNQ1(4) dNQ2(4)];
    dm12 = NLM.m12(q1i, q2i).*[dNQ1(1) dNQ2(1)] + ...
        NLM.m12(q1i, q2i+1)*[dNQ1(2) dNQ2(2)] + ...
        NLM.m12(q1i+1, q2i).*[dNQ1(3) dNQ2(3)] + ...
        NLM.m12(q1i+1, q2i+1)*[dNQ1(4) dNQ2(4)];
    
    dknl = dknl(:,1).*dQ1 + dknl(:,2).*dQ2;
    dcnl = dcnl(:,1).*dQ1 + dcnl(:,2).*dQ2;
    dm12 = dm12(1,1).*dQ1 + dm12(1,2).*dQ2;
    
    % forcing
    fm = [frc(1:2, q1i)*N1(1)+frc(1:2, q1i+1)*N1(2); 
        frc(3:4, q2i)*N2(1)+frc(3:4, q2i+1)*N2(2)];
    dfm = diag([frc(1:2, q1i)*dN1(1)+frc(1:2, q1i+1)*dN1(2); frc(3:4, q2i)*dN2(1)+frc(3:4, q2i+1)*dN2(2)]);
    
    dfm = dfm(:,1).*dQ1 + dfm(:,2).*dQ2;
    
    % Assembling
%     M = [1 m12; m12 1];
%     K = [knl(1) knl(2); knl(2) knl(3)];
%     C = [cnl(1) cnl(2); cnl(2) cnl(3)];
    
%     E = [K-w^2*M w*C;-w*C K-w^2*M];
    
%     R = [(K-w^2*M)*U(1:2) + (w*C)*U(3:4) - fm(1:2);
%         - (w*C)*U(1:2) + (K-w^2*M)*U(3:4) - fm(3:4)];
    
    R = [(knl(1)-w^2)    *U(1) + (knl(2)-w^2*m12)*U(2) + (w*cnl(1))      *U(3) + (w*cnl(2))      *U(4) - fm(1);
         (knl(2)-w^2*m12)*U(1) + (knl(3)-w^2)    *U(2) + (w*cnl(2))      *U(3) + (w*cnl(3))      *U(4) - fm(2);
        -(w*cnl(1))      *U(1) - (w*cnl(2))      *U(2) + (knl(1)-w^2)    *U(3) + (knl(2)-w^2*m12)*U(4) - fm(3);
        -(w*cnl(2))      *U(1) - (w*cnl(3))      *U(2) + (knl(2)-w^2*m12)*U(3) + (knl(3)-w^2)    *U(4) - fm(4)];
    
    dRdw = [(-2*w)    *U(1) + (-2*w*m12)*U(2) + (cnl(1))  *U(3) + (cnl(2))  *U(4) - fm(1);
            (-2*w*m12)*U(1) + (-2*w)    *U(2) + (cnl(2))  *U(3) + (cnl(3))  *U(4) - fm(2);
           -(cnl(1))  *U(1) - (cnl(2))  *U(2) + (-2*w)    *U(3) + (-2*w*m12)*U(4) - fm(3);
           -(cnl(2))  *U(1) - (cnl(3))  *U(2) + (-2*w*m12)*U(3) + (-2*w)    *U(4) - fm(4)];
       
    dRdU = [(dknl(1,:)-w^2)     *U(1) + (dknl(2,:)-w^2*dm12)*U(2) + (w*dcnl(1,:))       *U(3) + (w*dcnl(2,:))       *U(4) - dfm(1,:);
            (dknl(2,:)-w^2*dm12)*U(1) + (dknl(3,:)-w^2)     *U(2) + (w*dcnl(2,:))       *U(3) + (w*dcnl(3,:))       *U(4) - dfm(2,:);
           -(w*dcnl(1,:))       *U(1) - (w*dcnl(2,:))       *U(2) + (dknl(1,:)-w^2)     *U(3) + (dknl(2,:)-w^2*dm12)*U(4) - dfm(3,:);
           -(w*dcnl(2,:))       *U(1) - (w*dcnl(3,:))       *U(2) + (dknl(2,:)-w^2*dm12)*U(3) + (dknl(3,:)-w^2)     *U(4) - dfm(4,:)] + ...
           [(knl(1)-w^2)    , (knl(2)-w^2*m12), (w*cnl(1))      , (w*cnl(2));
            (knl(2)-w^2*m12), (knl(3)-w^2)    , (w*cnl(2))      , (w*cnl(3));
           -(w*cnl(1))      , -(w*cnl(2))     , (knl(1)-w^2)    , (knl(2)-w^2*m12);
           -(w*cnl(2))      , -(w*cnl(3))     , (knl(2)-w^2*m12), (knl(3)-w^2)];
    
    dRdUw = [dRdU dRdw];
end