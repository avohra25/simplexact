\ A staffing plan that cannot work. `simplexact iis examples/staffing.lp`
\ finds the few rows and bounds that conflict among the ones that don't matter.
minimize
  cost: 120 senior + 80 junior + 60 intern
subject to
  coverage: senior + junior + intern >= 12
  mentoring: 2 senior - intern >= 0
  budget: 120 senior + 80 junior + 60 intern <= 900
  review: senior + junior >= 6
  desks: senior + junior + intern <= 20
  juniors_need_seniors: junior - 3 senior <= 0
bounds
  senior <= 3
  junior <= 8
  intern <= 10
end
