\ Beale's example: plain Dantzig pivoting cycles forever on this degenerate LP.
\ The default rule falls back to Bland's rule and reaches the optimum -5/4.
minimize
  -3/4 a + 20 b - 1/2 c + 6 d
subject to
  r1: 1/4 a - 8 b - c + 9 d <= 0
  r2: 1/2 a - 12 b - 1/2 c + 3 d <= 0
  r3: c <= 1
end
