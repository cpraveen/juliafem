module fem

using Data
using PyPlot

function stima3(vertices)
   d = size(vertices,2);
   G = [ ones(1,d+1); vertices' ] \ [ zeros(1,d); eye(d) ];
   return det([ones(1,d+1); vertices']) * G * G' / prod(1:d);
end

function fem_50()
   coord = readdlm("coordinates.dat");
   elem3 = isfile("elements3.dat") ? round(Int64,readdlm("elements3.dat")) : [];
   elem4 = isfile("elements4.dat") ? round(Int64,readdlm("elements4.dat")) : [];
   neumann = isfile("neumann.dat") ? round(Int64,readdlm("neumann.dat")) : [];
   dirichlet = round(Int64,readdlm("dirichlet.dat"));
   n = size(coord,1); nt = size(elem3,1); nq = size(elem4,1);
   println("Number of vertices, triangles, quads = $n, $nt, $nq")

   A = sparse(Int64[],Int64[],Float64[],n,n);
   b = zeros(n);
   for j in 1:nt
      v = vec(elem3[j,:]);
      A[v,v] += stima3(coord[v,:]);
      b[v] += det([[1,1,1]'; coord[v,:]']) * f(sum(coord[v,:],1)/3)/6;
   end

   if ~isempty(neumann)
      nn = size(neumann,1);
      for j in 1:nn
         v = vec(neumann[j,:]);
         b[v] += norm(coord[v[1],:]-coord[v[2],:]) * g(sum(coord[v,:],1)/2)/2;
      end
   end

   u = zeros(n);
   BoundNodes = unique(dirichlet);
   u[BoundNodes] = u_d( coord[BoundNodes,:] );
   b -= A * u;

   FreeNodes = setdiff( 1:n, BoundNodes );
   u[FreeNodes] = A[FreeNodes,FreeNodes] \ b[FreeNodes];

   triplot(coord[:,1], coord[:,2], elem3-1)
   tricontour(coord[:,1], coord[:,2], elem3-1, u)
   show()
end

export fem_50

end
