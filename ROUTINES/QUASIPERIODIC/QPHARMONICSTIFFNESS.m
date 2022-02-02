function [E,dEdw] = QPHARMONICSTIFFNESS(M, C, K, w, h)
%QPHARMONICSTIFFNESS
%
%   USAGE :
%       M,C,K   : (Nd,Nd)
%       w       : (1,Nc)
%       h       : (Nh,Nc)
%   OUTPUT :
%       E       : (Nd*Nhc,Nd*Nhc)

    if sum(all(h==0, 2)) && ~all(h(1,:)==0)
        error('If you want the dc term, put it in the beginning of h')
    end

    Nc = size(h,2);
    Nh = size(h,1);
    Nhc = sum(all(h==0, 2)+2*any(h~=0, 2));

    hiw  = h*w(:);
    dhiw = h;
    D = kron(diag(hiw), [0 1;-1 0]);
    if all(h(1,:)==0)
        D = D(2:end, 2:end);
    end        
    
    E = kron(speye(Nhc),K) + kron(D^2,M) + kron(D, C);
    
    if nargout==2
        dEdw = cell(Nc,1);
        for ic=1:Nc
            dD = kron(diag(h(:,ic)), [0 1;-1 0]);
            if all(h(1,:)==0)
               dD = dD(2:end,2:end);
            end
            dEdw{ic} = -2*kron(D*dD, M) + kron(dD, C);
        end
        dEdw = cat(3, dEdw{:});
    end
end