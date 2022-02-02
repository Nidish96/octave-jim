function [y] = QPTIMETRANS(Y, h, Nt)
%QPTIMETRANS
%
%   USAGE :
%       Y       : (Nhc,Ny)
%       h       : (Nh,Nc)
%       Nt      :
%   OUTPUT :
%       y       : (Nt^Nc,Ny)

    if sum(all(h==0, 2)) && ~all(h(1,:)==0)
        error('If you want the dc term, put it in the beginning of h')
    end

    Nc = size(h,2);
    Nh = size(h,1);
    Nhc = sum(all(h==0, 2)+2*any(h~=0, 2));
    
    iof = fix(Nt/2)+1;
    i0 = num2cell(repmat(iof, 1, Nc));
   
    Ny = size(Y, 2);
    y = zeros(Nt^Nc,Ny);
    for yi=1:Ny
        yf = zeros(repmat(Nt, 1, Nc));
        k = 1;
        for hi=1:Nh
            if all(h(hi,:)==0)
                yf(i0{:}) = Y(k,yi)*2;
                k = k+1;
            else
                id = num2cell(iof+h(hi,:));
                yf(id{:}) = Y(k,yi)-1j*Y(k+1,yi);
                id = num2cell(iof-h(hi,:));
                yf(id{:}) = Y(k,yi)+1j*Y(k+1,yi);
                k = k+2;
            end
        end
        yf = yf*(Nt^Nc)/2;

        % IFFT
        y(:, yi) = reshape(ifftn(ifftshift(yf)), Nt^Nc, 1);
    end
end