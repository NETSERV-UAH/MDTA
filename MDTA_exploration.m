
function [distance2source] = MDTA_exploration(netCostMatrix, s)
%==============================================================
% shortestPath: the list of nodes in the shortestPath from source to destination;
% totalCost: the total cost of the  shortestPath;
% farthestNode: the farthest node to reach for each node after performing the routing;
% n: the number of nodes in the network;
% s: source node index;
% d: destination node index;
%==============================================================
%  Original Dijsktra's algorithm Code by:
% ++by Xiaodong Wang
% ++23 Jul 2004 (Updated 29 Jul 2004)
% ++http://www.mathworks.com/matlabcentral/fileexchange/5550-dijkstra-shortest-path-routing
% Modifications (simplifications) by Meral Shirazipour 9 Dec 2009
%==============================================================
%
% Modified by Diego López Pajares, Universidad de Alcalá
for j=1:size(netCostMatrix,3)
    netCostMatrix_aux=netCostMatrix(:,:,j);
    n = size(netCostMatrix_aux,1);
    visited(1:n) = false;
    distance(1:n) = inf;
    distance2source_aux= Inf(size(netCostMatrix_aux));
    parent(1:n) = 0;
    distance(s) = 0;
    for i = 1:(n)
        temp = [];
        for h = 1:n
            if ~visited(h)
                temp=[temp distance(h)];
            else
                temp=[temp inf];
            end
        end
        [t, u] = min(temp);
        visited(u) = true;
        for v = 1:n
            if(v~=s)
                if(parent(u) ~= v)
                    if ( ( netCostMatrix_aux(u, v) + distance(u)) < distance(v) )
                        distance(v) = distance(u) + netCostMatrix_aux(u, v);
                        parent(v) = u;
                    end
                    distance2source_aux(u,v)= distance(u) + netCostMatrix_aux(u, v);
                else
                    distance2source_aux(u,v)= distance(u) + netCostMatrix_aux(u,v);
                end
            end
        end
    end
    distance2source(:,:,j)=distance2source_aux;
end
end