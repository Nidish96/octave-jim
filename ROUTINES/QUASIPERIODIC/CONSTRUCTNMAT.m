function [Nmat, varargout] = CONSTRUCTNMAT(ws, Nc, Nt, varargin)  % Need to think about w-gradients
%CONSTRUCTNMAT constructs the "Nmat" interpolation matrix for the different
%QP marching approaches
%
%   USAGE:
%       [Nmat] = CONSTRUCTNMAT(ws, Nc, Nt, type);  % for type=1,2,4
%           (OR)
%       [Nmat, bpis, bpjs] = CONSTRUCTNMAT(ws, Nc, Nt, 3);  % for type=3

    if nargin==3
        Nmtype = 1;
    else
        Nmtype = varargin{1};
    end
    deltau = 2*pi/Nt; 
    taus = cell(Nc, 1);
    [taus{:}] = ndgrid(1:Nt);    

    tausnm1 = cell(Nc-1, 1);
    [tausnm1{:}] = ndgrid(1:Nt);    
    switch Nmtype
        case 1
            dt_vec = ws*deltau/vecnorm(ws);  % vector corresponding to deltatau amplitude in real time dxn
            dt_vec_dw = permute(deltau*(vecnorm(ws)*eye(Nc)-ws(:).*ws/vecnorm(ws))/vecnorm(ws)^2, [3 2 1]);  % derivative of dt_vec wrt ws

            ptsb = cellfun(@(c) str2double(c), num2cell(dec2bin(0:(2^Nc-1))));  % using binary for the construction of points on a unit square
            oppi = (2^Nc):-1:1;    % diagonally opposite points are retrieved using the binary inverses
            xis = ptsb*(-deltau);  % coordinates in "tau" space relative to origin

            Lm = deltau^Nc;                   % Lebesgue Measure of each cell in tau space
            rcoords = xis(oppi,:)-(-dt_vec);  % relative coordinates from the "opposite points" on each cell
            Nsf = prod(abs(rcoords), 2)'/Lm;  % Lagrange Shape Functions to interpolate previous point in the cell diagonally behind (in tau space)

            % Compute derivate of Nsf wrt ws
            rcoords_dw = (-(-dt_vec_dw)).*sign(rcoords);
            Nsf_dw = permute(sum(cell2mat(arrayfun(@(a) (prod(abs(rcoords(:,setdiff(1:Nc,a))),2).*rcoords_dw(:, a, :))/Lm, 1:Nc, 'UniformOutput', false)), 2), [2 1 3]);

            ijs = fix(cell2mat(cellfun(@(c) c(:), taus, 'UniformOutput', false)'));  % indices of all points
            evijs = mod(repmat(ijs, 1, 1, 2^Nc) - permute(ptsb, [3 2 1])-1, Nt)+1;   % indices of points forming the cell that is diagonally behind
            evns = squeeze(sum((Nt.^(0:Nc-1)).*(evijs(:, 1:end, :)-1),2)+1);         % vertex IDs (according to list in ijs) corresponding to points forming the cell that is diagonally behind

            % Build sparse interpolation matrix
            Nmat = sparse(repmat((1:Nt^Nc)', 1, 2^Nc), evns, repmat(Nsf, Nt^Nc, 1));

            % Sort Nodes by distance from origin
            ds = mod(plus(taus{:})-1-1,Nt)+1;
            bpis = cell(Nt,1);
            bpjs = cell(Nt,1);
            for ti=1:Nt
                bpis{ti} = find(ds==ti);
            end
            for ti=1:Nt
                tim1 = mod(ti-1-1, Nt)+1;
                tim2 = mod(ti-2-1, Nt)+1;
                bpjs{ti} = unique([bpis{tim1}; bpis{tim2}]);  % Unique not necessary.
            end
            varargout{1} = bpis;
            varargout{2} = bpjs;
            if nargout==4  % Also return omega-gradient of Nmat if necessary
                varargout{3} = arrayfun(@(a) sparse(repmat((1:Nt^Nc)', 1, 2^Nc), evns, repmat(Nsf_dw(:, :, a), Nt^Nc, 1)), 1:Nc, 'UniformOutput', false);
            end
        case 2
            ptsb = eye(Nc);

            Nsf = ws/sum(ws);  % Shape functions constructed using omegas
            Nsf_dw = permute(eye(Nc)/sum(ws)-ws/sum(ws)^2, [3 2 1]);

            ijs = cell2mat(cellfun(@(c) c(:), taus, 'UniformOutput', false)');  % indices of all points
            evijs = mod(repmat(ijs, 1, 1, Nc) - permute(ptsb, [3 2 1])-1, Nt)+1;  % indices of points forming the FD stencil
            evns = squeeze(sum((Nt.^(0:Nc-1)).*(evijs(:, 1:end, :)-1),2)+1);  % vertex IDs (according to list in ijs) corresponding to points from the stencil

            % Build sparse interpolation matrix
            Nmat = sparse(repmat((1:Nt^Nc)', 1, Nc), evns, repmat(Nsf, Nt^Nc, 1));

            % Sort Nodes by distance from origin
            ds = mod(plus(taus{:})-1-1,Nt)+1;
            bpis = cell(Nt,1);
            for ti=1:Nt
                bpis{ti} = find(ds==ti);
            end
            bpjs = bpis([Nt 1:Nt-1]);
            varargout{1} = bpis;
            varargout{2} = bpjs;
            if nargout==4  % Also return omega-gradient of Nmat if necessary
                varargout{3} = arrayfun(@(a) sparse(repmat((1:Nt^Nc)', 1, Nc), evns, repmat(Nsf_dw(:, :, a), Nt^Nc, 1)), 1:Nc, 'UniformOutput', false);
            end
        case 3  % Problems Problems. Nmat(bpis{1}, bpis{1}) has a few non-zero elements.
            ptsb = ones(1, Nc);  % diagonally opposite point
            
            ijs = cell2mat(cellfun(@(c) c(:), taus, 'UniformOutput', false)');  % indices of all points
            evijs = mod(repmat(ijs, 1, 1, 1) - permute(ptsb, [3 2 1])-1, Nt)+1;  % indices of points forming the cell that is diagonally behind
            
            evns = squeeze(sum((Nt.^(0:Nc-1)).*(evijs(:, 1:end, :)-1),2)+1);  % vertex IDs (according to list in ijs) corresponding to points forming the cell that is diagonally behind
            
            % Build interpolant matrix
            dt_vecs = ws * deltau./ws(:);  % Each row is the "previous point" projected on a different plane. We will have to apply the correct one for the correct points
            dt_vecs_dw = deltau*( ws(:).*permute(eye(Nc),[3 1 2])-ws.*permute(eye(Nc), [1 3 2]) )./ws(:).^2;  % derivatives of the above 
            
            ptsbm = cellfun(@(c) str2double(c), num2cell(dec2bin(0:(2^(Nc-1)-1))));  % using binary for the construction of points on a unit square
            Lm = deltau^(Nc-1);
            
            inds = reshape((1:Nt^Nc)', [repmat(Nt, 1, Nc) ones(Nc==1)]);
            S.subs = repmat({':'}, 1, Nc);
            S.type = '()';
            bpis = cell(Nt,1);  bpjs = cell(Nt,1);    vals = cell(Nt,1);
            vals_dw = cell(Nt,1);
            for ti=1:Nt  % march over the diagonal
                bpis{ti} = zeros((Nt-ti+1)^(Nc-1)*Nc,1);
                bpjs{ti} = zeros((Nt-ti+1)^(Nc-1)*Nc,2^(Nc-1));
                vals{ti} = zeros((Nt-ti+1)^(Nc-1)*Nc,2^(Nc-1));
                vals_dw{ti} = zeros((Nt-ti+1)^(Nc-1)*Nc,2^(Nc-1), Nc);
                for ci=1:Nc
                    S.subs = repmat({ti:Nt}, 1, Nc);
                    S.subs{ci} = ti;
                    sinds = unique(reshape(subsref(inds, S), [], 1));
                
                    bpis{ti}((ci-1)*(Nt-ti+1)^(Nc-1)+(1:(Nt-ti+1)^(Nc-1))) = sinds;
            
                    % build local Lagrange interpolant for these points
                    cinds = setxor(1:Nc, ci);
                    dt_loc = -dt_vecs(ci, cinds);
                    dt_loc_dw = -dt_vecs_dw(ci, cinds, :);
            
                    n1s_loc = floor(dt_loc/deltau);  % relevant point on the diagonal
                    
                    xis = (n1s_loc+ptsbm)*deltau;  % points on a cube around previous projection
            
                    rcoords = xis(end:-1:1,:)-dt_loc;  % relative coordinates
                    Nsfs = prod(abs(rcoords),2)'/Lm;  % Lagrangian Shape Functions
                    
                    rcoords_dw = -dt_loc_dw.*sign(rcoords);
                    if Nc>2
                        Nsfs_dw = permute(sum(cell2mat(arrayfun(@(a) rcoords_dw(:, a, :).*rcoords(:, setdiff(1:Nc-1, a)), 1:Nc-1, 'UniformOutput', false)), 2), [2 1 3]);
                    elseif Nc==2
                        Nsfs_dw = permute(rcoords_dw, [2 1 3])/Lm;
                    end

                    ptsb = zeros(2^(Nc-1), Nc);
                    ptsb(:, cinds) = n1s_loc+ptsbm+1;
            
                    bpevijs = mod(repmat(ijs(evns(sinds),:), 1, 1, 2^(Nc-1)) + permute(ptsb, [3 2 1])-1, Nt)+1;  % indices of points forming the relevant cell
                    bpjs{ti}((ci-1)*(Nt-ti+1)^(Nc-1)+(1:(Nt-ti+1)^(Nc-1)),:) = squeeze(sum((Nt.^(0:Nc-1)).*(bpevijs(:, 1:end, :)-1),2)+1);  % vertex IDs (according to list in ijs) 
                    vals{ti}((ci-1)*(Nt-ti+1)^(Nc-1)+(1:(Nt-ti+1)^(Nc-1)),:) = repmat(Nsfs, (Nt-ti+1)^(Nc-1), 1);

                    vals_dw{ti}((ci-1)*(Nt-ti+1)^(Nc-1)+(1:(Nt-ti+1)^(Nc-1)),:,:) = reshape(cell2mat(arrayfun(@(a) repmat(Nsfs_dw(:, :, a), (Nt-ti+1)^(Nc-1), 1), 1:Nc, 'UniformOutput', false)), (Nt-ti+1)^(Nc-1), 2^(Nc-1), Nc);
                end
                [bpis{ti}, indtis] = unique(bpis{ti});
                bpjs{ti} = bpjs{ti}(indtis,:);
                vals{ti} = vals{ti}(indtis,:);
                vals_dw{ti} = vals_dw{ti}(indtis,:,:);
            end
            % Correction for the corner points on the first set
            [~, si] = intersect(bpis{1}, bpjs{1});
            qi = [];  qj = [];
            for i=1:length(si)
                [i, j] = find(bpjs{1}==bpis{1}(si(i)));
                qi = [qi;i];
                qj = [qj;j];
            end
            idealbpj1 = setdiff(bpjs{1}, bpis{1});
            for i=1:length(qi)  % Fix interpolants of problematic points 
                dis = bpjs{1}(qi(i),:);  % indices of nodes. qj(i)^th node has to be replaced.
                d_coords = cell2mat(cellfun(@(t) t(dis), taus, 'UniformOutput', false))';  % Each row is a node
                pp_coord = vals{1}(qi(i),:)*d_coords;
                pp_coord_dw = reshape(cell2mat(arrayfun(@(a) vals_dw{1}(qi(i),:,a)*d_coords, 1:Nc, 'UniformOutput',false)), 1, Nc, Nc);

                dns = [];
                si = [];
                its = 1;
                while isempty(dns) || isempty(si)
                    dns = mod(dis(qj(i))+kron(its*[-1 1],Nt.^(0:Nc-1))-1, Nt^Nc)+1;
                    % Eliminate points already considered in interpolation
                    ois = setdiff(1:Nc,qj(i));
                    [~, si] = intersect(dns, dis(ois));
                    dns(si) = [];  % Already considered in current interpolation
                    [~, si] = intersect(dns, bpis{1});
                    dns(si) = [];
                    [~, si] = intersect(dns, idealbpj1);
                    its = its+1;
                end
                dns = dns(si(1));  % Choose just another point

                dis(qj(i)) = dns;
                dn_coords = cell2mat(cellfun(@(t) t(dis), taus, 'UniformOutput', false))';  % Each row is a node
                dcp = d_coords(qj(i),:)-pp_coord;
                dcpn = dn_coords(qj(i),:)-pp_coord;
                for ic=find(dcp~=0)
                    if sign(dcpn(ic))~=sign(dcp(ic))
                        dn_coords(qj(i), ic) = dn_coords(qj(i), ic) + sign(dcp(ic))*(Nt);
                    end
                end
                wgts = pp_coord/dn_coords;
                wgts_dw = reshape(cell2mat(arrayfun(@(a) pp_coord_dw(:, :, a)/dn_coords, 1:Nc, 'UniformOutput', false)), 1, 2^(Nc-1), Nc);

                bpjs{1}(qi(i),:) = dis;
                vals{1}(qi(i),:) = wgts;
                vals_dw{1}(qi(i),:,:) = wgts_dw;
            end

            varargout{1} = bpis;
            varargout{2} = bpjs;
            [bpis, uinds] = unique(cell2mat(bpis));
            bpjs = cell2mat(bpjs);  bpjs = bpjs(uinds,:);
            vals = cell2mat(vals);  vals = vals(uinds,:);
            vals_dw = cell2mat(vals_dw); vals_dw = reshape(vals_dw(uinds,:), Nt^Nc, 2^(Nc-1), Nc);
            
            % Build sparse interpolation matrix
            Nmat = sparse(repmat(bpis, 1, 2^(Nc-1)), bpjs, vals);
            if nargout==4  % frequency gradient of Nmat if necessary
                varargout{3} = arrayfun(@(a) sparse(repmat(bpis, 1, 2^(Nc-1)), bpjs, vals_dw(:, :, a)), 1:Nc, 'UniformOutput', false);
            end
        case 4   % March on Stagger-Stretched Grid
            dt_vec = ws*deltau/ws(end);
            dt_vec_dw = permute(eye(Nc), [3 1 2])*deltau/ws(end);
            dt_vec_dw(:, :, end) = deltau*([zeros(1,Nc-1) 1]/ws(end)-ws/ws(end)^2);
            
            lp_in = floor(dt_vec(1:end-1)*Nt/deltau)+1;  % "Lower corner" on initial plane for interpolation
            
            ptsb = cellfun(@(c) str2double(c), num2cell(dec2bin(0:(2^(Nc-1)-1))));  % using binary for the construction of points on a unit square
            oppi = (2^(Nc-1)):-1:1;  % diagonally opposite points are retrieved using the binary inverses
            
            xis = ptsb*deltau;  % Cell coordinates relative to origin point
            
            ptc = dt_vec(1:end-1)*Nt-(lp_in-1)*deltau';  % relative coordinate of where origin from initial plane will wrap around on the initial plane
            ptc_dw = dt_vec_dw(1,1:end-1,:)*Nt;  % w-gradient of above
            
            Lm = deltau^(Nc-1);  % Lebesgue Measure of each cell in tau space
            rcoords = xis(oppi,:)-ptc;
            Nsf = prod(abs(rcoords), 2)'/Lm;  % Lagrange Shape Functions to interpolate previous point in the cell diagonally behind (in tau space)

            rcoords_dw = -ptc_dw.*sign(rcoords);
            if Nc>2
                Nsf_dw = permute(sum(cell2mat(arrayfun(@(a) rcoords_dw(:, a, :).*rcoords(:, setdiff(1:Nc-1, a)), 1:Nc-1, 'UniformOutput', false)), 2), [2 1 3])/Lm;
            elseif Nc==2
                Nsf_dw = permute(rcoords_dw, [2 1 3])/Lm;
            end
            
            ijs = cell2mat(cellfun(@(c) c(:), tausnm1, 'UniformOutput', false)');  % List of indices on (Nc-1)D grid signifying points on "next plane"
            
            fp_ijs = mod((ijs-1)+lp_in-1, Nt)+1;  % List of lower vertex of containing cell on "initial plane"
            evijs = mod(repmat(fp_ijs, 1, 1, 2^(Nc-1)) + permute(ptsb, [3 2 1])-1, Nt)+1;
            evns = squeeze(sum((Nt.^(0:Nc-2)).*(evijs(:, 1:end, :)-1),2)+1);
            
            % UNLIKE the Nmat in the other types, this only interpolates
            % the "next plane" using the "initial plane"
            Nmat = sparse( repmat((1:Nt^(Nc-1))', 1, 2^(Nc-1)), evns, repmat(Nsf, Nt^(Nc-1), 1) );
            if nargout==2
                varargout{1} = arrayfun(@(a) sparse( repmat((1:Nt^(Nc-1))', 1, 2^(Nc-1)), evns, repmat(Nsf_dw(:, :, a), Nt^(Nc-1), 1) ), 1:Nc, 'UniformOutput', false);
            end
        case 5  % Same as case 4 but we compute the inverse mapping (inverse of matrix Nmat) directly. Not working unfortunately, idk why :(3 jan 22)
            dt_vec = ws*deltau/ws(end);
            dt_vec_dw = permute(eye(Nc), [3 1 2])*deltau/ws(end);
            dt_vec_dw(:, :, end) = deltau*([zeros(1,Nc-1) 1]/ws(end)-ws/ws(end)^2);
                        
            lp_ne = floor(-dt_vec(1:end-1)*Nt/deltau);
            
            ptsb = cellfun(@(c) str2double(c), num2cell(dec2bin(0:(2^(Nc-1)-1))));  % using binary for the construction of points on a unit square
            oppi = (2^(Nc-1)):-1:1;  % diagonally opposite points are retrieved using the binary inverses
            
            xis = ptsb*deltau;  % Coordinates relative to origin point
            
            ptc = -dt_vec(1:end-1)*Nt-(lp_ne)*deltau;  % relative coordinate of "next point" in cell
            ptc_dw = -dt_vec_dw(1,1:end-1,:)*Nt;  % w-gradient of above
            
            Lm = deltau^(Nc-1);  % Lebesgue Measure of each cell in tau space
            rcoords = xis(oppi,:)-ptc;
            Nsf = prod(abs(rcoords), 2)'/Lm;  % Lagrange Shape Functions to interpolate previous point in the cell diagonally behind (in tau space)

            rcoords_dw = -ptc_dw.*sign(rcoords);
            if Nc>2
                Nsf_dw = permute(sum(cell2mat(arrayfun(@(a) rcoords_dw(:, a, :).*rcoords(:, setdiff(1:Nc-1, a)), 1:Nc-1, 'UniformOutput', false)), 2), [2 1 3])/Lm;
            elseif Nc==2
                Nsf_dw = permute(rcoords_dw, [2 1 3])/Lm;
            end
            
            ijs = cell2mat(cellfun(@(c) c(:), tausnm1, 'UniformOutput', false)');  % List of indices on (Nc-1)D grid signifying points on "initial plane"
            
            ne_ijs = mod((ijs-1)+lp_ne, Nt)+1;  % List of lower vertex of containing cell on "next plane"
            evijs = mod(repmat(ne_ijs, 1, 1, 2^(Nc-1)) + permute(ptsb, [3 2 1])-1, Nt)+1;
            evns = squeeze(sum((Nt.^(0:Nc-2)).*(evijs(:, 1:end, :)-1),2)+1);
         
            % UNLIKE the Nmat in the other types, this only interpolates
            % the "initial plane" using the "next plane"
            Nmat = sparse( repmat((1:Nt^(Nc-1))', 1, 2^(Nc-1)), evns, repmat(Nsf, Nt^(Nc-1), 1) );
            if nargout==2
                varargout{1} = arrayfun(@(a) sparse( repmat((1:Nt^(Nc-1))', 1, 2^(Nc-1)), evns, repmat(Nsf_dw(:, :, a), Nt^(Nc-1), 1) ), 1:Nc, 'UniformOutput', false);
            end
    end
end