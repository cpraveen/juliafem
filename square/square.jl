using PyPlot

n = 10;
println("n = $n");

xmin=0.0; xmax=1.0;
ymin=0.0; ymax=1.0;
nx=n;
ny=n;
x=LinRange(xmin,xmax,nx);
y=LinRange(ymin,ymax,ny);
c=1;
X=zeros(nx*ny,1);
Y=zeros(nx*ny,1);
for j in 1:nx
   for k in 1:ny
      X[c] = x[j];
      Y[c] = y[k];
      global c = c + 1;
   end
end

tri=DelaunayTri([X,Y]);
triplot(tri)
e=tri.edges;           # all edges
eb=tri.freeBoundary;   # boundary edges

np   = length(X);   # no of vertices
nt   = size(tri,1); # no of triangles
println("Number of points         = $np");
println("Number of triangles      = $nt");

n_e  = size(e ,1); # no of all edges
n_eb = size(eb,1); # no of boundary edges

# find neighbouring cell of boundary edges
for j in 1:n_eb
   nbr = tri.edgeAttachments(eb(j,:));
   nbr = nbr{:};
   dx = X(eb(j,2)) - X(eb(j,1));
   dy = Y(eb(j,2)) - Y(eb(j,1));
   # order vertices so that domain is to left side
   v  = tri.Triangulation(nbr(1),:);
   dx1= sum(X(v))/3 - X(eb(j,1));
   dy1= sum(Y(v))/3 - Y(eb(j,1));
   if(dx1*dy - dx*dy1 > 0)
      tmp = eb(j,:);
      eb(j,1) = tmp(2);
      eb(j,2) = tmp(1);
   end
end

writedlm("coordinates.dat",tri.X');
writedlm("elements3.dat", tri.Triangulation');
writedlm("dirichlet.dat", eb);

#fid=fopen('neumann.dat','w');
#fprintf(fid,'%8d %8d\n', eb);
#fclose(fid);
