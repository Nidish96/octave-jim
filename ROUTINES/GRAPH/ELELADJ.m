function [Ac, Am] = ELELADJ(MESH)
%NODEELADJ Creates an adjacency matrix describing the element
%connections in the given 4-noded 2D mesh. The weights will be
%the distances from the centroids.
% USAGE:
%	[A] = NODEELADJ(MESH);
% INPUTS:
%   MESH 	: MESH Structure with,
% 	Nds, Tri, Quad, Nn, Ne
% OUTPUTS:
%   Ac		: NexNe adjacency matrix with centroidal distances as
%   		  weights
%   Am		: NexNe adjacency matrix with type of connection as 
% 		  weights

    Nn = MESH.Nn;
    Ne = MESH.Ne;

    Ac = zeros(Ne, Ne);
    Am = zeros(Ne, Ne);
    
    [Q1, ~] = ZTE_ND2QP(MESH, 1);
    Ctrds = Q1*MESH.Nds;  % Element Centroids
    
    % Elements Sharing 1 node: distance = 2
    % Elements Sharing 2 nodes: distance = 1
    for e=1:Ne
        for i=1:4
            nf = find(sum(MESH.Quad(:, 2:end)==MESH.Quad(e, i+1),2));
        
            Ac(e, nf) = Ac(e, nf) + 1;
        end
        Ac(e, e) = 0;
        nf = find(Ac(e,:));
        Am(e, nf) = vecnorm(Ctrds(nf,:)-Ctrds(e,:), 2, 2);
    end
end