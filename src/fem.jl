module fem

using Data
using PyPlot

function stima3(vertices)
   d = size(vertices,2)
   G = [ ones(1,d+1); vertices' ] \ [ zeros(1,d); eye(d) ]
   return det([ones(1,d+1); vertices']) * G * G' / prod(1:d)
end

function stima4(vertices)
   D_Phi = [ vertices[2,:]-vertices[1,:] vertices[4,:]-vertices[1,:] ]
   B = inv( D_Phi' * D_Phi )
   C1 =  [ [ 2, -2] [-2,  2] ] * B[1,1]
       + [ [ 3,  0] [ 0, -3] ] * B[1,2]
       + [ [ 2,  1] [ 1,  2] ] * B[2,2]
   C2 =  [ [-1,  1] [ 1, -1] ] * B[1,1]
       + [ [-3,  0] [ 0,  3] ] * B[1,2]
       + [ [-1, -2] [-2, -1] ] * B[2,2]
  return det( D_Phi ) * [ C1 C2; C2 C1 ] / 6
end

function plotsol(coord, elem3, elem4, u, ncont)
   triplot(coord[:,1], coord[:,2], elem3-1, color=(0.0,0.25,0.15),linewidth=0.2)
   tricontour(coord[:,1], coord[:,2], elem3-1, u, ncont, linewidth=2)

   # For quads, we triangulate and plot
   if size(elem4,1) > 0
      el = elem4[:,[1,2,3]]-1
      triplot(coord[:,1], coord[:,2], el, color=(0.0,0.25,0.15),linewidth=0.2)
      tricontour(coord[:,1], coord[:,2], el, u, ncont, linewidth=2)
      el = elem4[:,[1,3,4]]-1
      triplot(coord[:,1], coord[:,2], el, color=(0.0,0.25,0.15),linewidth=0.2)
      tricontour(coord[:,1], coord[:,2], el, u, ncont, linewidth=2)
   end
   show()
end

function fem_50(ncont=10)
   coord = readdlm("coordinates.dat")
   elem3 = isfile("elements3.dat") ? round(Int64,readdlm("elements3.dat")) : []
   elem4 = isfile("elements4.dat") ? round(Int64,readdlm("elements4.dat")) : []
   neumann = isfile("neumann.dat") ? round(Int64,readdlm("neumann.dat")) : []
   dirichlet = round(Int64,readdlm("dirichlet.dat"))
   n, nt, nq = size(coord,1), size(elem3,1), size(elem4,1)
   println("Number of vertices, triangles, quads = $n, $nt, $nq")

   A = spzeros(n,n)
   b = zeros(n)

   for j in 1:nt
      v = vec(elem3[j,:])
      A[v,v] += stima3(coord[v,:])
      b[v] += det([[1,1,1]'; coord[v,:]']) * f(sum(coord[v,:],1)/3)/6
   end

   for j in 1:nq
      v = vec(elem4[j,:])
      A[v,v] += stima4(coord[v,:])
      b[v] += det([[1,1,1]'; coord[v[1:3],:]']) * f(sum(coord[v,:],1)/4)/4
   end

   if ~isempty(neumann)
      nn = size(neumann,1)
      for j in 1:nn
         v = vec(neumann[j,:])
         b[v] += norm(coord[v[1],:]-coord[v[2],:]) * g(sum(coord[v,:],1)/2)/2
      end
   end

   u = zeros(n)
   BoundNodes = unique(dirichlet)
   u[BoundNodes] = u_d( coord[BoundNodes,:] )
   b -= A * u

   FreeNodes = setdiff( 1:n, BoundNodes )
   u[FreeNodes] = A[FreeNodes,FreeNodes] \ b[FreeNodes]

   plotsol(coord, elem3, elem4, u, ncont)
end

export fem_50

end
