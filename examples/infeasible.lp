\ No point satisfies all three rows; the solver explains why with Farkas multipliers.
maximize
  x + y
subject to
  a: x + 2 y <= 4
  b: 3 x + y <= 6
  c: x + y >= 5
end
