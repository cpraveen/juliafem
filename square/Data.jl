module Data

# This returns a scalar
function g(u)
   return 0.0;
end

function u_d(u)
  return zeros( size( u, 1 ), 1 );
end

# This returns a scalar
function f(u)
   return 2.0*pi*pi*sin(pi*u[1])*sin(pi*u[2]);
end

export g, u_d, f

end
