module Data

# Neumann bc, this returns a scalar
function g(u)
   return 0.0
end

# Dirichlet bc
function u_d(u)
  return zeros( size( u, 1 ), 1 )
end

# RHS function in poisson equation. This returns a scalar
function f(u)
   return 1.0
end

export g, u_d, f

end
