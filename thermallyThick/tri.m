
% Simple tridiagonal solver using Thomas Algorithm
% --                         -- --      --   --     --
% | b_1 c_1                   | | fout_1 | = | fin_1 |
% | a_2 b_2 c_2               | | fout_2 | = | fin_2 |
% |      a_3 b_3 c_3          | | fout_3 | = | fin_3 |
% |        .   .   .          | |    .   | = |   .   |
% |            .   .    .     | |    .   | = |   .   |
% |                .    .   . | |    .   | = |   .   |
% |                   a_M b_M | | fout_M | = | fin_M |
% --                         -- --      --   --     --
% Note: x is a dummy array
%        f array is overwritten by solution
%        a, b and c are preserved
%
function f = tri(a,b,c,f)
  M = length(a);
  x = zeros(M,1);
  x(1) = c(1)/b(1);
  f(1) = f(1)/b(1);
  % forward sweep
  for ji=2:M
     z = 1./( b(ji) - a(ji)*x(ji-1) );
     x(ji) = c(ji)*z;
     f(ji) = ( f(ji) - a(ji)*f(ji-1) )*z;
  end
  % backward sweep
  for ji=M-1:-1:1
     f(ji) = f(ji)- x(ji)*f(ji+1);
  end
end
