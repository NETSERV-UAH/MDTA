function [path_list,netCostMatrix_updated,visited_nodes,trees]=MDTA_construction(netCostMatrix, src, mode,dst, n_paths,varargin)
% Created by Diego López Pajares, Universidad de Alcalá
% MDTA_construction returns multiple disjoint paths between the source
% node and one or multiple destination nodes, considering the information
% that already has the net cost matrix, which has been analyzed using
% Dijkstra's algorithm. The matrix could be two or three-dimensional,
% depending on whether it is the first time calculation.
% Each two-dimensional plane represents the differt topology analysis for each path.
% The disjointness is achieved by labeling the paths used in one plane
% with an infinite cost in the remaining planes. This labeling prevents other
% paths from using it.
% input:
%     - netCostMatrix: cost matrix. Dim1 y Dim2: main matrix; Dim3: disjoint
%       path plane. first time iteration has only two dimensions with the original
%       cost matrix in each iterations new planes will be mounted over the main
%       one
%     - src: source node
%     - mode: disjoint mode "node" disjoint, "link" disjoint
%     - dst: destiantion nodes array
%     - n_paths: max number of paths that that will be attempted
% Output
%     - path_list: discovered paths list
%     - netCostMatrix_updated: cost matrix updated with the discovered paths
netCostMatrix_updated=netCostMatrix;
n_dst=length(dst);
path_list=cell(1,n_dst);
% Operation mode -> identifier
% 1 -> link
% 2 -> node
if(strcmp(mode,'link'))
    op_mode=1;
elseif (strcmp(mode,'node'))
    op_mode=2;
elseif (strcmp(mode,'nodeIterative'))
    op_mode=3;
else
    error('Operation mode possiblities are link or node, review input paramenters');
end
if(nargin==6)
    visited_nodes=varargin{1};
else
    visited_nodes=[];
end
for i=1:n_dst
    switch(op_mode)
        case 1
            [path_list_dst,netCostMatrix_updated]=LinkMode(netCostMatrix_updated,src,dst(i),n_paths);
            if isempty(path_list_dst)
                path_list(1,i)={[]};
            else
                path_list(1,i)=path_list_dst;
            end
        case 2
            netCostMatrix_orig=netCostMatrix_updated;
            if((ismember(dst(i),visited_nodes)))
                netCostMatrix_updated=MDP_RestoreMtrxCostForNewDstNode(dst(i),...
                    netCostMatrix,netCostMatrix_updated,visited_nodes);
            end
            netCostMatrix_copy=netCostMatrix_updated;
            [path_list_dst,netCostMatrix_updated,visited_nodes]=...
                NodeMode(netCostMatrix_updated,src,dst(i),n_paths,visited_nodes);
            netCostMatrix_updated=MDP_RestoreMtrxCostOnlyNewLinksDiscovered(dst(i),visited_nodes,path_list_dst,netCostMatrix_orig,netCostMatrix_copy,netCostMatrix_updated);
            if isempty(path_list_dst)
                path_list(1,i)={[]};
            else
                path_list(1,i)=path_list_dst;
            end
    end
end
trees=cell(n_paths,1);
path_list_ref=[path_list{:,:}];
counters = cell2mat(cellfun(@(cell) cellfun(@numel, cell), path_list, 'UniformOutput', false));
counters=counters-(counters~=0);
adders=(counters~=0);
idx=adders;
while(any(counters~=0,"all"))
    slide=cellfun(@slideWind,path_list_ref,num2cell(double(idx)),'UniformOutput',false);
    slide=cell2mat(reshape(slide',[],1));
    slide = mat2cell(slide, repmat(size(slide, 1)/n_paths, 1, n_paths), size(slide, 2));
    trees = cellfun(@treeupdate,slide,trees,'UniformOutput',false);
    counters=counters-adders;
    adders=(counters~=0);
    idx=(idx+adders).*adders;
end
end
function [path_list, new_mtrx_cost]=LinkMode(mtrx_cost,src,dst,n_paths)
path_list=cell(1,1);
n_paths_heritaged=size(mtrx_cost,3);
new_path_search=false;
prev_hop=0;
if(n_paths_heritaged==1)
    n_paths_heritaged=n_paths;
    new_mtrx_cost=ones(size(mtrx_cost,1),size(mtrx_cost,2),n_paths_heritaged);
    new_mtrx_cost=new_mtrx_cost.*mtrx_cost;
    new_path_search=true;
else
    new_mtrx_cost=mtrx_cost;
end

for i=1:n_paths_heritaged
    mtrx_cost_aux=new_mtrx_cost;
    [sorted,nodes]=sort(mtrx_cost_aux(:,:,i),1,'ascend');
    nodes(sorted==inf)=nan;
    nodes_copy=nodes;
    prev_hop=0;
    path=[dst];
    next_hop=dst;
    while (next_hop~=src && ~isnan(next_hop))
        [prev_hop,next_hop,nodes_copy]=GetNextNode(next_hop,nodes_copy);
        if(prev_hop==dst && isnan(next_hop))
            break;
        end
        if((next_hop==dst || isnan(next_hop) || ismember(next_hop,path)))
            [prev_hop,next_hop,path]=CheckHopBack(prev_hop, path, nodes_copy);
        else
            [mtrx_cost_aux]= UpdateMtrxCost(next_hop,prev_hop, mtrx_cost_aux,i);
            path=[path next_hop];
        end
        if(next_hop==src)

            path=flip(path);
            [mtrx_cost_aux,path]=CorrectLoopInPath(path,mtrx_cost,mtrx_cost_aux);
            new_mtrx_cost=mtrx_cost_aux;
            path_list{1,1}{end+1,1}=(path);
        end
    end
    if(next_hop~=src)
        path_list{1,1}{end+1,1}=[];
    end
end

end
function [prv_h, nxt_h, node_indx]= GetNextNode(current_node,nodes_list)
prv_h=current_node;
nxt_h=nodes_list(1,current_node);
nodes_list(:,prv_h)=[nodes_list(2:end,prv_h); nan];
if(~isnan(nxt_h))
    [row,cols,val]=find(nodes_list(:,nxt_h)==prv_h);
    nodes_list(row:end,nxt_h)=[nodes_list(row+1:end,nxt_h);nan];
end
node_indx=nodes_list;
end
function [prv_h, nxt_h, path]=CheckHopBack(prv_hop, path, nodes_list)
duplicate=find(path==prv_hop);
if(isnan(nodes_list(1,path(duplicate(end)))))
    duplicate=duplicate-1;
end
if(duplicate(end)==0)
    prv_h=path(1);
    path=path(1);
else
    prv_h=path(duplicate(end));
    path=path(1:duplicate(end));
end
nxt_h=prv_h;
end

function [mtrx]= UpdateMtrxCost(next_hop,prv_host, mtrx, iter)
mtrx(next_hop,prv_host,iter)=mtrx(next_hop,prv_host,iter)/1000000; %Reduce cost reverse
mtrx(prv_host,next_hop,iter)=mtrx(prv_host,next_hop,iter)/10000; %Reduce cost
mtrx(next_hop,prv_host,setdiff(1:end,iter))=inf;
mtrx(prv_host,next_hop,setdiff(1:end,iter))=inf;
end
function [mtrxCostModified,new_path]=CorrectLoopInPath(path,originalMtrxCost,mtrxCostModified)
[~,order1]=unique(path,'first');
[~,order2]=unique(path,'last');
loops_detected=[setdiff(order1,order2,'stable'), setdiff(order2,order1,'stable')];
if(~isempty(loops_detected))
    length_of_the_loop=loops_detected(:,2)-loops_detected(:,1);
    [~,dim_max_loop]=max(length_of_the_loop);
    ranges2remove=~((loops_detected(:,1)>loops_detected(dim_max_loop,1) & loops_detected(:,1)<loops_detected(dim_max_loop,2)) | (loops_detected(:,2)>loops_detected(dim_max_loop,1) & loops_detected(:,2)<loops_detected(dim_max_loop,2)));
    loops_2_delete=loops_detected(ranges2remove,:);
    loops_2_delete(:,2)=loops_2_delete(:,2)+1;
    loops_2_delete=sortrows(loops_2_delete,2);
    nodes_deleted=cell(size(loops_2_delete,1));
    for i=1:size(loops_2_delete,1)
        if(i==1)
            new_path=path(1:loops_2_delete(i,1));
        else
            new_path=[new_path,path(loops_2_delete(i-1,2):loops_2_delete(i,1))];
        end
        nodes_deleted{i}=path(loops_2_delete(i,1):loops_2_delete(i,2)-1);
    end
    new_path=[new_path,path(loops_2_delete(i,2):end)];
    for i=1:length(nodes_deleted)
        for j=1:length(nodes_deleted{i})-1
            x=nodes_deleted{i}(j);
            y=nodes_deleted{i}(j+1);
            mtrxCostModified(x,y,:)=originalMtrxCost(x,y,:);
            mtrxCostModified(y,x,:)=originalMtrxCost(y,x,:);
        end
    end
else
    new_path=path;
end
end

function [path_list, new_mtrx_cost,visited_nodes,visited_nodes_end]=NodeMode(mtrx_cost,src,dst,n_paths,visited_nodes)
path_list=cell(1,1);
new_path_search=false;
n_paths_heritaged=size(mtrx_cost,3);
if(n_paths_heritaged==1)
    n_paths_heritaged=n_paths;
    new_mtrx_cost=ones(size(mtrx_cost,1),size(mtrx_cost,2),n_paths_heritaged);
    new_mtrx_cost=new_mtrx_cost.*mtrx_cost;
    new_path_search=true;
else
    new_mtrx_cost=mtrx_cost;
end

for i=1:n_paths_heritaged
    prev_hop=0;
    mtrx_cost_aux=new_mtrx_cost;
    [sorted,nodes]=sort(mtrx_cost_aux(:,:,i),1,'ascend');
    nodes(sorted==inf)=nan;
    nodes_copy=nodes;
    path=[dst];
    next_hop=dst;
    while (next_hop~=src && ~isnan(next_hop))
        [prev_hop,next_hop,nodes_copy]=GetNextNodeNodeMode(next_hop,nodes_copy,path);
        if(prev_hop==dst && isnan(next_hop))
            break;
        end
        if((next_hop==dst || isnan(next_hop)))
            prev_hop_before=prev_hop;
            [prev_hop,next_hop,path,mtrx_cost_aux]=CheckHopBackNodeMode(prev_hop, path, nodes_copy,mtrx_cost_aux,new_mtrx_cost);
            visited_nodes(visited_nodes==prev_hop_before)=[];

        else
            [mtrx_cost_aux]= UpdateMtrxCostNodeMode(next_hop,prev_hop, mtrx_cost_aux,i,src); % with no hop back update costs
            path=[path next_hop];
            if(~ismember(next_hop,visited_nodes))
                visited_nodes=[visited_nodes next_hop];
            end
        end
        if(next_hop==src)
            path=flip(path);
            new_mtrx_cost=mtrx_cost_aux;
            path_list{1,1}{end+1,1}=(path);
        end
    end
    if(next_hop~=src)
        path_list{1,1}{end+1,1}=[];
    end
end
visited_nodes_end=visited_nodes;
end
function [prv_h, nxt_h, node_indx]= GetNextNodeNodeMode(current_node,nodes_list,path)
prv_h=current_node;
i=1;
nxt_h=nodes_list(i,current_node);
while(ismember(nxt_h,path) && ~isnan(nxt_h))
    i=i+1;
    nxt_h=nodes_list(i,current_node);
end
nodes_list(:,prv_h)=[nodes_list(2:end,prv_h); nan];
if(~isnan(nxt_h))
    [row,cols,val]=find(nodes_list(:,nxt_h)==prv_h);
    nodes_list(row:end,nxt_h)=[nodes_list(row+1:end,nxt_h);nan];
end
node_indx=nodes_list;
end

function [mtrx]= UpdateMtrxCostNodeMode(next_hop,prv_host, mtrx, iter,src)
mtrx(next_hop,prv_host,iter)=mtrx(next_hop,prv_host,iter)/1000000;
mtrx(prv_host,next_hop,iter)=mtrx(prv_host,next_hop,iter)/10000;
if(next_hop~=src)
    mtrx(next_hop,:,setdiff(1:end,iter))=inf;
    mtrx(:,next_hop,setdiff(1:end,iter))=inf;
    mtrx(next_hop,prv_host,setdiff(1:end,iter))=inf;
    mtrx(prv_host,next_hop,setdiff(1:end,iter))=inf;
end
end
function [prv_h, nxt_h, path,mtrx]=CheckHopBackNodeMode(prv_hop, path, nodes_list,mtrx,prev_mtrx)
duplicate=find(path==prv_hop);
if(isnan(nodes_list(1,path(duplicate(end)))))
    duplicate=duplicate-1;
end
if(duplicate(end)==0)
    prv_h=path(1);
    path=path(1);
else
    prv_h=path(duplicate(end));
    path=path(1:duplicate(end));
end
nxt_h=prv_h;
mtrx(prv_h,:,:)=prev_mtrx(prv_h,:,:);
mtrx(:,prv_h,:)=prev_mtrx(:,prv_h,:);
end
function outp=slideWind(inp,idx)
outp=[0,0];
if(idx~=0)
    outp=inp(idx:idx+1);
end
end
function tree=treeupdate(slide,tree)
slide=unique(slide,"rows");
if(isempty(tree))
    tree=zeros(0,2) ;
end
if(any(~ismember(slide,tree,"rows") & ~ all(slide==0,2)))
    tree=[tree;slide(~ismember(slide,tree,"rows") & ~ all(slide==0,2),:)];
end

end