\ Stigler-style diet problem: meet nutrient minimums at the lowest cost.
\ Costs are in cents per serving, so the optimum comes out as an exact fraction.
minimize
  cost: 18 oats + 23 milk + 5 bread + 41 eggs + 30 beans
subject to
  calories: 110 oats + 160 milk + 70 bread + 150 eggs + 120 beans >= 2000
  protein:  4 oats + 8 milk + 3 bread + 12 eggs + 8 beans >= 55
  calcium:  2 oats + 285 milk + 22 bread + 54 eggs + 80 beans >= 800
  iron:     2/5 oats + 1/10 milk + 3/5 bread + 9/5 eggs + 5/2 beans >= 15
bounds
  oats <= 4
  milk <= 8
  bread <= 12
  eggs <= 2
  beans <= 4
end
