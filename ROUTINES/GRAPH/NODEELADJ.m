function [Ad, An] = NODEELADJ(Nds,Els)
%NODEELADJ Creates an adjacency matrix describing the nodal
%connections in the given 4-noded 2D FEM model. The weights will be
%the nodal distances.
% USAGE:
%	[Ad] = NODEELADJ(Nds,Els);
% INPUTS:
%   Nds		: Ndx2 matrix of nodal coordinates
%   Els		: Nex4 matrix describing elements (and nodes
%   		  contained therein)
% OUTPUTS:
%   Ad		: NdxNd adjacency matrix with nodal distances as
%   		  weights
%   An		: NdxNd adjacency matrix with number of connections
%                 as weights

    Nd = size(Nds,1);
    Ne = size(Els,1);
    
    Ad = zeros(Nd,Nd);
    An = zeros(Nd,Nd);
    
    V = zeros(4,2);
    Vd = zeros(4,1);
    for e=1:Ne
        V = Nds(Els(e,:),:);
        Vd = sqrt(sum((diff(V([1:end 1],:))).^2,2));
        
        Ad(Els(e,1),Els(e,2)) = Vd(1);
        Ad(Els(e,2),Els(e,3)) = Vd(2);
        Ad(Els(e,3),Els(e,4)) = Vd(3);
        Ad(Els(e,4),Els(e,1)) = Vd(4);
        
        Ad(Els(e,2),Els(e,1)) = Vd(1);
        Ad(Els(e,3),Els(e,2)) = Vd(2);
        Ad(Els(e,4),Els(e,3)) = Vd(3);
        Ad(Els(e,1),Els(e,4)) = Vd(4);
        
        for i=1:4
            An(Els(e,i),Els(e,i)) = An(Els(e,i),Els(e,i)) + 1;
        end
        
        An(Els(e,1),Els(e,2)) = An(Els(e,1),Els(e,2)) + 1;
        An(Els(e,2),Els(e,3)) = An(Els(e,2),Els(e,3)) + 1;
        An(Els(e,3),Els(e,4)) = An(Els(e,3),Els(e,4)) + 1;
        An(Els(e,4),Els(e,1)) = An(Els(e,4),Els(e,1)) + 1;
        
        An(Els(e,2),Els(e,1)) = An(Els(e,2),Els(e,1)) + 1;
        An(Els(e,3),Els(e,2)) = An(Els(e,3),Els(e,2)) + 1;
        An(Els(e,4),Els(e,3)) = An(Els(e,4),Els(e,3)) + 1;
        An(Els(e,1),Els(e,4)) = An(Els(e,1),Els(e,4)) + 1;
    end    
end