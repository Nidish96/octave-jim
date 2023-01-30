function [Ns, varargout] = QPTIMEINTERP(taus, h, varargin)
%QPTIMEINTERP
%
%   USAGE :
%       taus    : (Np, Nc)
%       h       : (Nh, Nc)
%   OUTPUTS :
%       Ns      : (Np, Nhc)
    if sum(all(h==0, 2)) && ~all(h(1,:)==0)
        error('If you want the dc term, put it in the beginning of h')
    end

    Nc = size(h,2);
    Nh = size(h,1);
    Nhc = sum(all(h==0, 2)+2*any(h~=0, 2));
    Np = size(taus, 1);
    
    Ns = zeros(Np, Nhc);
    k = 1;
    if all(h(1,:)==0)
        Ns(:, k) = 1;
        k = k+1;
    end
    Ns(:, k:2:end) = cos(taus*h(k:end,:)');
    Ns(:, k+1:2:end) = sin(taus*h(k:end,:)');

    if nargout==2
        taus_dw = varargin{1};
        Ns_dw = zeros(Np, Nhc, Nc);
        k = 1 + all(h(1,:)==0);
        Ns_dw(:, k:2:end, :) = -sin(taus*h(k:end,:)').*cell2mat(permute(arrayfun(@(a) taus_dw(:,:,a)*h(k:end,:)',1:Nc,'UniformOutput', false), [1 3 2]));
        Ns_dw(:, k+1:2:end, :) = cos(taus*h(k:end,:)').*cell2mat(permute(arrayfun(@(a) taus_dw(:,:,a)*h(k:end,:)',1:Nc,'UniformOutput', false), [1 3 2]));

        varargout{1} = Ns_dw;
    end
end