module Data

# Neumann bc, this returns a scalar
function g(u)
   return 0.0;
end

# Dirichlet bc
function u_d(u)
  return zeros( size( u, 1 ), 1 );
end

# RHS function in poisson equation. This returns a scalar
function f(u)
   return 2.0*pi*pi*sin(pi*u[1])*sin(pi*u[2]);
end

export g, u_d, f

end
