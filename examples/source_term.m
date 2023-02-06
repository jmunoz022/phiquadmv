function ff0 = source_term(source,t0,M,nel,nnode,coord,nodes,nquad,dir)

%%%Source vector%%%
ff = zeros(nnode,1);
for iel = 1:nel
    nd = nodes(iel,:);
    xcoord = coord(nd,1);
    ycoord = coord(nd,2);
    
    f = local_source(source,nquad,xcoord,ycoord,t0);
    for i = 1:length(nd)
        ii = nd(i);
        ff(ii) = ff(ii)+f(i);
    end
end 

%%%homogeneous Dirchlet BC%%%
ff(dir) = [];
ff0 = M\ff;


end