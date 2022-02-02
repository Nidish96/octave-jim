function [Y] = QPFOURIERCOEFF_old(y, h)
%QPFOURIERCOEFF
%
%   USAGE :
%       y       : (Nt,Nt,...) or (Ny,1) cell of them
%       h       : (Nh,Nc)
%   OUTPUT :
%       Y       : (Nhc,Ny)

    if sum(all(h==0, 2)) && ~all(h(1,:)==0)
        error('If you want the dc term, put it in the beginning of h')
    end

    Nc = size(h,2);
    Nh = size(h,1);
    Nhc = sum(all(h==0, 2)+2*any(h~=0, 2));
    
    if ~strcmp(class(y),'class')
        y = {y};
    end
    Ny = length(y);
    Nt = size(y{1}, 1);
    
    iof = fix(Nt/2)+1;
    i0 = num2cell(repmat(iof, 1, Nc));
    
    Y = zeros(Nhc, Ny);
    for yi=1:Ny
        % Fourier Transform & Scale
        yf = fftshift(fftn(y{yi}))*2/(Nt^Nc);
        yf(i0{:}) = yf(i0{:})/2;

        % Choose selected Harmonics
        k = 1;
        for hi=1:Nh
            if all(h(hi,:)==0)  % zero harmonics
                Y(k,yi) = yf(i0{:});
                k = k+1;
            else
                id = num2cell(iof+h(hi,:));
                Y(k,yi)   = real(yf(id{:}));
                Y(k+1,yi) = -imag(yf(id{:}));
                k = k+2;
            end
        end
    end
end