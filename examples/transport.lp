\ Transportation problem: two plants ship to three markets at minimum cost.
minimize
  cost: 2.5 sea_ny + 1.7 sea_chi + 1.8 sea_top + 2.5 sd_ny + 1.8 sd_chi + 1.4 sd_top
subject to
  sea_supply: sea_ny + sea_chi + sea_top <= 350
  sd_supply:  sd_ny + sd_chi + sd_top <= 600
  ny_demand:  sea_ny + sd_ny >= 325
  chi_demand: sea_chi + sd_chi >= 300
  top_demand: sea_top + sd_top >= 275
end
