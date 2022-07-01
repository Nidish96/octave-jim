function [Nmat, varargout] = CONSTRUCTNMATHERE(ws, Nc, Nt, varargin)  % Need to think about w-gradients
%CONSTRUCTNMATHERE constructs the Nmat
    if nargin==3
        Nmtype = 1;
    else
        Nmtype = varargin{1};
    end
    deltau = 2*pi/Nt; 
    taus = cell(Nc, 1);
    [taus{:}] = ndgrid(1:Nt);    
    switch Nmtype
        case 1
            dt_vec = ws*deltau/vecnorm(ws);  % vector corresponding to deltatau amplitude in real time dxn

            ptsb = cellfun(@(c) str2double(c), num2cell(dec2bin(0:(2^Nc-1))));  % using binary for the construction of points on a unit square
            oppi = (2^Nc):-1:1;    % diagonally opposite points are retrieved using the binary inverses
            xis = ptsb*(-deltau);  % coordinates in "tau" space relative to origin

            Lm = deltau^Nc;                                 % Lebesgue Measure of each cell in tau space
            Nsf = prod(abs(xis(oppi,:)-(-dt_vec)), 2)'/Lm;  % Lagrange Shape Functions to interpolate previous point in the cell diagonally behind (in tau space)

            ijs = fix(cell2mat(cellfun(@(c) c(:), taus, 'UniformOutput', false)'));  % indices of all points
            evijs = mod(repmat(ijs, 1, 1, 2^Nc) - permute(ptsb, [3 2 1])-1, Nt)+1;   % indices of points forming the cell that is diagonally behind
            evns = squeeze(sum((Nt.^(0:Nc-1)).*(evijs(:, 1:end, :)-1),2)+1);         % vertex IDs (according to list in ijs) corresponding to points forming the cell that is diagonally behind

            % Build sparse interpolation matrix
            Nmat = sparse(repmat((1:Nt^Nc)', 1, 2^Nc), evns, repmat(Nsf, Nt^Nc, 1));
        case 2
            ptsb = eye(Nc);

            Nsf = ws/sum(ws);  % Shape functions constructed using omegas

            ijs = cell2mat(cellfun(@(c) c(:), taus, 'UniformOutput', false)');  % indices of all points
            evijs = mod(repmat(ijs, 1, 1, Nc) - permute(ptsb, [3 2 1])-1, Nt)+1;  % indices of points forming the FD stencil
            evns = squeeze(sum((Nt.^(0:Nc-1)).*(evijs(:, 1:end, :)-1),2)+1);  % vertex IDs (according to list in ijs) corresponding to points from the stencil

            % Build sparse interpolation matrix
            Nmat = sparse(repmat((1:Nt^Nc)', 1, Nc), evns, repmat(Nsf, Nt^Nc, 1));
        case 3
            ptsb = ones(1, Nc);  % diagonally opposite point
            
            ijs = cell2mat(cellfun(@(c) c(:), taus, 'UniformOutput', false)');  % indices of all points
            evijs = mod(repmat(ijs, 1, 1, 1) - permute(ptsb, [3 2 1])-1, Nt)+1;  % indices of points forming the cell that is diagonally behind
            
            evns = squeeze(sum((Nt.^(0:Nc-1)).*(evijs(:, 1:end, :)-1),2)+1);  % vertex IDs (according to list in ijs) corresponding to points forming the cell that is diagonally behind
            
            % Build interpolant matrix
            dt_vecs = ws * deltau./ws(:);  % Each row is the "previous point" projected on a different plane. We will have to apply the correct one for the correct points
            
            ptsbm = cellfun(@(c) str2double(c), num2cell(dec2bin(0:(2^(Nc-1)-1))));  % using binary for the construction of points on a unit square
            Lm = deltau^(Nc-1);
            
            inds = reshape((1:Nt^Nc)', [repmat(Nt, 1, Nc) ones(Nc==1)]);
            S.subs = repmat({':'}, 1, Nc);
            S.type = '()';
            bpis = cell(Nt,1);  bpjs = cell(Nt,1);    vals = cell(Nt,1);
            for ti=1:Nt  % march over the diagonal
                bpis{ti} = zeros((Nt-ti+1)^(Nc-1)*Nc,1);
                bpjs{ti} = zeros((Nt-ti+1)^(Nc-1)*Nc,2^(Nc-1));
                vals{ti} = zeros((Nt-ti+1)^(Nc-1)*Nc,2^(Nc-1));
                for ci=1:Nc
                    S.subs = repmat({ti:Nt}, 1, Nc);
                    S.subs{ci} = ti;
                    sinds = unique(reshape(subsref(inds, S), [], 1));
                
                    bpis{ti}((ci-1)*(Nt-ti+1)^(Nc-1)+(1:(Nt-ti+1)^(Nc-1))) = sinds;
            
                    % build local Lagrange interpolant for these points
                    cinds = setxor(1:Nc, ci);
                    dt_loc = -dt_vecs(ci, cinds);
            
                    n1s_loc = floor(dt_loc/deltau);  % relevant point on the diagonal
                    
                    xis = (n1s_loc+ptsbm)*deltau;  % points on a cube around previous projection
            
                    Nsfs = abs(prod(xis(end:-1:1,:)-dt_loc,2)'/Lm);  % Lagrangian Shape Functions
                    
                    ptsb = zeros(2^(Nc-1), Nc);
                    ptsb(:, cinds) = n1s_loc+ptsbm+1;
            
                    bpevijs = mod(repmat(ijs(evns(sinds),:), 1, 1, 2^(Nc-1)) + permute(ptsb, [3 2 1])-1, Nt)+1;  % indices of points forming the relevant cell
                    bpjs{ti}((ci-1)*(Nt-ti+1)^(Nc-1)+(1:(Nt-ti+1)^(Nc-1)),:) = squeeze(sum((Nt.^(0:Nc-1)).*(bpevijs(:, 1:end, :)-1),2)+1);  % vertex IDs (according to list in ijs) 
                    vals{ti}((ci-1)*(Nt-ti+1)^(Nc-1)+(1:(Nt-ti+1)^(Nc-1)),:) = repmat(Nsfs, (Nt-ti+1)^(Nc-1), 1);
                end
                [bpis{ti}, indtis] = unique(bpis{ti});
                bpjs{ti} = bpjs{ti}(indtis,:);
                vals{ti} = vals{ti}(indtis,:);
            end
            varargout{1} = bpis;
            varargout{2} = bpjs;
            [bpis, uinds] = unique(cell2mat(bpis));
            bpjs = cell2mat(bpjs);  bpjs = bpjs(uinds,:);
            vals = cell2mat(vals);  vals = vals(uinds,:);
            
            % Build sparse interpolation matrix
            Nmat = sparse(repmat(bpis, 1, 2^(Nc-1)), bpjs, vals);
    end
end