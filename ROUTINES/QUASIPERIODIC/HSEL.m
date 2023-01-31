function [h] = HSEL(Nhmax, ws, varargin)
    if isempty(varargin)
        hcr = 1;
    else
        hcr = varargin{1};
    end
    Nc = length(ws);

    %% Harmonic Selection
    hall = cell(1, Nc);
    [hall{:}] = ndgrid(-Nhmax:Nhmax);
    hall = cell2mat(cellfun(@(c) c(:), hall, 'UniformOutput', false));
    
    switch hcr
        case 1 % Criterion 1  - Abs Criterion
            h = hall(sum(abs(hall),2)<=Nhmax & sum(hall,2)>=0,:);
            
            % Legacy
            h(sum(h,2)==0 & h(:,1)<=0, :) = [];
            h = [zeros(1,Nc); h];
            % Remove when obsolete
        case 2 % Criterion 2 - w-Criterion
            hpos = hall;
            hws = hpos*ws(:);
            [~, si] = sort(hws);
            h = hpos(si(1:Nhmax^Nc),:);
        otherwise
            error('hcr should be one of [1,2]. %d given.', hcr)
    end
%     h(sum(h,2)==0 & h(:,1)<=0, :) = [];
%     for ic=1:Nc
%         id = find(all(h == [zeros(1,ic-1) 1 zeros(1, Nc-ic)], 2));
%         if ~isempty(id)
%             h(id,:) = [];
%         end
%     end
%     h = [zeros(1,Nc); eye(Nc); h];
end