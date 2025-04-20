using JLD2
@load "Desktop/JMMenura/paper_scripts/example_scripts/parameter_selection/anole_sim_diff/AlphaOUAnoles.jld2" alphas wts
wtsOU = wts

@load "Desktop/JMMenura/paper_scripts/example_scripts/parameter_selection/anole_sim_diff/AlphaOUISOAnoles.jld2" pars wts
wtsISO = wts

@load "/home/simoneb/Desktop/JMMenura/paper_scripts/example_scripts/parameter_selection/anole_sim_diff/BManoles.jld2" sigmas wts

wtsBM = wts

sum(sum(wtsISO))/(sum(sum(wtsOU)) + sum(sum(wtsBM)) + sum(sum(wtsISO)))  # 0.6344563749450002
sum(sum(wtsBM))/(sum(sum(wtsOU)) + sum(sum(wtsBM)) + sum(sum(wtsISO)))  # 0.32664209482954354
sum(sum(wtsOU))/(sum(sum(wtsOU)) + sum(sum(wtsBM)) + sum(sum(wtsISO))) # 0.03890153022545631