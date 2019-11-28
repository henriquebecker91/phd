library(dplyr)
library(reshape2)
library(xtable)

get_experiments <- function () {
  final <- read.table("./csv/experiment0.csv", header=T, sep=";", stringsAsFactors = FALSE)
  final$experiment <- 0
  for (e in c(1, 2, 3, 17)) {
    fpath = paste0("./csv/experiment", e, ".csv")
    tmp <- read.table(fpath, header=T, sep=";", stringsAsFactors = FALSE)
    tmp$experiment <- e
    final <- bind_rows(final, tmp)
  }
  for (e in c(4, 5, 6, 7, 8, 18)) {
    fpath = paste0("./csv/experiment_", e, ".csv")
    tmp <- read.table(fpath, header=T, sep=";", stringsAsFactors = FALSE)
    tmp$experiment <- e
    final <- bind_rows(final, tmp)
  }
  # experiments 9 to 12
  tmp <- read.table("./csv/martin_furini.csv", header=T, sep=";", stringsAsFactors = FALSE)
  final <- bind_rows(final, tmp)
  # experiments 13 to 16
  tmp <- read.table("./csv/dimitri_furini_abs.csv", header=T, sep=";", stringsAsFactors = FALSE)
  final <- bind_rows(final, tmp)
  # experiment 20
  tmp <- read.table("./csv/dimitri_furini_priced_abs.csv", header=T, sep=";", stringsAsFactors = FALSE)
  final <- bind_rows(final, tmp)
  
  final$time_to_solve_model[final$stop_code != 1] <- NA
  final <- final %>% arrange(experiment, instfname, seed)
  final$instfname <- sapply(final$instfname, basename)
  
  return(final)
}

bignumber2latex <- function (x) {
  expo <- trunc(log10(x))
  nume <- round(x / (10^expo), digits = 3)
  cat(paste0("\\(", nume, " \\times 10^", expo, "\\)"))
}

t <- get_experiments()

# experiments description:
# 0: just model building, faithful, no reductions
# 1: just model building, faithful, only cut position
# 2: just model building, faithful, only redundant cut
# 3: just model building, faithful, both reductions
# 4: just seven easy, faithful, no reductions
# 5: just seven easy, faithful, both reductions
# 6: just seven easy, faithful, both reductions and round2disc
# 7: all 33 + 5, revised, both reductions* 
# 8: all 33 + 5, revised, both reductions* and round2disc
# 9: just model building, martin reimplementation, no reductions
# 10: just model building, martin reimplementation, only cut position
# 11: just model building, martin reimplementation, only redundant cut
# 12: just model building, martin reimplementation, both reductions
# 13: just extracted data, original implementation, no reductions
# 14: just extracted data, original implementation, only cut position
# 15: just extracted data, original implementation, only redundant cut
# 16: just extracted data, original implementation, both reductions
# 17: (RAN BUT NOT EXTRACTED) just model building, faithful, both reductions and round2disc
# 18: (RAN BUT NOT EXTRACTED) all 33 + 5, revised, no reduction
# 19: (ABANDONED IDEA) all 33 + 5, just model building, revised, no reductions but round2disc?
# 20: just extracted data, original implementation, number of variables in final priced version
# * Redundant-Cut is superseded, no?

# just checking if everything is correct:
# num_runs <- t %>% group_by(experiment) %>% tally()

# Analysis presented in text only:
# * how our re-implementation compares to the original in terms of plates and variables
#   + the number of plates is the same, the variables are not distant more than x%
t0 <- t %>% filter(experiment == 0 & !grepl("^P1_.*", instfname)) # our reimpl. of complete
t13 <- t %>% filter(experiment == 13) # extracted data from original article
stopifnot(length(t0$instfname) == length(t13$instfname))
stopifnot(sum(t0$num_plate == t13$num_plates) == 33)
stopifnot((sum(t0$num_vars)/sum(t13$num_vars)) < 1.03)
stopifnot(max(t0$num_vars/t13$num_vars) < 1.34)
stopifnot(min(t0$num_vars/t13$num_vars) < 1.01)
#t9 <- t %>% filter(experiment == 9) # extracted data from martin implementation
#t9$num_vars/t13$num_vars
# * how much we did reduce with cut position and redundant-cut, how much Martin reduce, how much furini reduced
t3 <- t %>% filter(experiment == 3 & !grepl("^P1_.*", instfname)) # furini reimplm. and both reductions
vars_0_3_per_red <- t3$num_vars/t0$num_vars
var_0_3_per_red <- round((1 - sum(t3$num_vars)/sum(t0$num_vars))*100, digits = 2)
plates_0_3_per_red <- round((1 - sum(t3$num_plates)/sum(t0$num_plates))*100, digits = 2)
#t12 <- t %>% filter(experiment == 12) 
#vars_9_12_per_red <- t12$num_vars/t9$num_vars
#abs(mean(vars_0_3_per_red) - mean(vars_9_12_per_red))*100 # less than 1% difference in reduction
#max(abs(vars_0_3_per_red - vars_9_12_per_red))
t16 <- t %>% filter(experiment == 16) 
vars_13_16_per_red <- t16$num_vars/t13$num_vars
var_13_16_per_red <- sum(t16$num_vars)/sum(t13$num_vars)
var_3_16_ratio <- sum(t3$num_vars)/sum(t16$num_vars)

var_13_16_per_red <- round((1 - sum(t16$num_vars)/sum(t13$num_vars))*100, digits = 2)
plates_13_16_per_red <- round((1 - sum(t16$num_plates)/sum(t13$num_plates))*100, digits = 2)

#abs(mean(vars_0_3_per_red) - mean(vars_13_16_per_red))*100
max(abs(vars_0_3_per_red - vars_13_16_per_red)*100)

# number of variables and plates of no-reduction revised model against original data
# with and without reductions
t18_33i_s1 <- t %>% filter(experiment == 18 & !grepl("^P1_.*", instfname) & seed == 1)
round((sum(t18_33i_s1$num_vars)/sum(t13$num_vars))*100, digits = 2)
sum(t18_33i_s1$num_vars)
round((sum(t18_33i_s1$num_plates)/sum(t13$num_plates))*100, digits = 2)
sum(t18_33i_s1$num_plates)
sum(t13$num_vars)
sum(t13$num_plates)
round((sum(t18_33i_s1$num_vars)/sum(t16$num_vars))*100, digits = 2)
round((sum(t18_33i_s1$num_plates)/sum(t16$num_plates))*100, digits = 2)
sum(t16$num_vars)
sum(t16$num_plates)
t7_33i_s1 <- t %>% filter(experiment == 7 & !grepl("^P1_.*", instfname) & seed == 1)
(1 - (sum(t7_33i_s1$num_vars)/sum(t18_33i_s1$num_vars)))*100
(1 - (sum(t7_33i_s1$num_plates)/sum(t18_33i_s1$num_plates)))*100
t8_33i_s1 <- t %>% filter(experiment == 8 & !grepl("^P1_.*", instfname) & seed == 1)
sum(t8_33i_s1$num_vars)
sum(t8_33i_s1$num_plates)
round((sum(t8_33i_s1$num_vars)/sum(t7_33i_s1$num_vars))*100, digits = 2)
round((sum(t8_33i_s1$num_plates)/sum(t7_33i_s1$num_plates))*100, digits = 2)

t7_33i <- t %>% filter(experiment == 7 & !grepl("^P1_.*", instfname))
t8_33i <- t %>% filter(experiment == 8 & !grepl("^P1_.*", instfname))
inst_solved_by_t7 <- t7_33i %>% group_by(instfname) %>% summarise(max_stop_code = max(stop_code))
inst_solved_by_t8 <- t8_33i %>% group_by(instfname) %>% summarise(max_stop_code = max(stop_code))

inst_solved_by_t7_t8 <- full_join(t7_33i, t8_33i, by = c("instfname", "seed")) %>%
  group_by(instfname) %>% summarise(max_stop_code = max(stop_code.x, stop_code.y))
inst_solved_by_t7_t8 %>% filter(max_stop_code > 1)
inst_solved_by_t7_t8 <- inst_solved_by_t7_t8 %>% filter(max_stop_code == 1)
inst_solved_by_t7_t8 <- inst_solved_by_t7_t8$instfname
t7_33i_solved <- t %>% filter(experiment == 7 & instfname %in% inst_solved_by_t7_t8)
t8_33i_solved <- t %>% filter(experiment == 8 & instfname %in% inst_solved_by_t7_t8)
t7_33i_solved_summ <- t7_33i_solved %>% group_by(instfname) %>% summarize(mean_time = mean(time_to_build_model + time_to_solve_model), sd_time = )
t8_33i_solved_summ <- t8_33i_solved %>% group_by(instfname) %>% summarize(mean_time = mean(time_to_build_model + time_to_solve_model))
sum(t7_33i_solved_summ$mean_time)
sum(t8_33i_solved_summ$mean_time)
sum(t8_33i_solved_summ$mean_time)/sum(t7_33i_solved_summ$mean_time)

t20 <- t %>% filter(experiment == 20)
round((sum(t8_33i_s1$num_vars)/sum(t20$num_vars))*100, digits = 2)

t17_33i <- t %>% filter(experiment == 17 & !grepl("^P1_.*", instfname))
sum(t3$num_vars)
sum(t17_33i$num_vars)
round((sum(t17_33i$num_vars)/sum(t3$num_vars))*100, digits = 2)
sum(t3$num_plates)
sum(t17_33i$num_plates)
round((sum(t17_33i$num_plates)/sum(t3$num_plates))*100, digits = 2)

# let us check the same for just the seven used in the experiments
t4 <- t %>% filter(experiment == 4) # get a experiment with just the seven instances
inst7 <- unique(t4$instfname)
t0_7i <- t0[t0$instfname %in% inst7, ]
t13_7i <- t13[t13$instfname %in% inst7, ]
t3_7i <- t3[t3$instfname %in% inst7, ]
t16_7i <- t16[t16$instfname %in% inst7, ]
sum(t0_7i$num_vars)/sum(t13_7i$num_vars)
sum(t3_7i$num_vars)/sum(t16_7i$num_vars)

# Tables
# * small table with five velasco instances
#   * lines: header + 5 instances
#   * columns: #instfname, #num_vars, #num_plates, #solve time, #lb, #ub
t8_5i_table <- t %>% filter(experiment == 8 & grepl("^P1_.*", instfname)) %>% group_by(instfname) %>%
  summarise(
    num_vars = dplyr::first(num_vars),
    num_plates = dplyr::first(num_plates),
    solve_time = mean(time_to_build_model + time_to_solve_model, na.rm = T),
    lb = max(obj_value),
    ub = min(obj_bound),
    finished = sum(stop_code == 1)
  )
t8_5i_table[is.nan(t8_5i_table$solve_time), ]$solve_time <- 3600
t8_5i_xtable <- xtable(t8_5i_table)
digits(t8_5i_xtable) <- c(0, 0, 0, 0, 2, 0, 0, 0)
print(t8_5i_xtable, only.contents = T, include.rownames = F, include.colnames = F)

# * small table at paper body comparing faithful both reductions against ours with both reduction
#   + lines: seven instances
#   + columns: (#vars, #plates, solve time, sd solve time, finished) * 2
#   + say that the build time was not more than x%? present total time?
#   + maybe say sd was always low?
t6 <- t %>% filter(experiment == 6)
inst7 <- unique(t6$instfname)
t8_7i <- t8_33i %>% filter(instfname %in% inst7)
t8_7i_summ <- t8_7i %>% group_by(instfname) %>% summarize(
  num_vars = dplyr::first(num_vars),
  num_plates = dplyr::first(num_plates),
  mean_time = mean(time_to_build_model + time_to_solve_model),
  sd_time = sd(time_to_build_model + time_to_solve_model)
)
t6_summ <- t6 %>% group_by(instfname) %>% summarize(
  num_vars = dplyr::first(num_vars),
  num_plates = dplyr::first(num_plates),
  mean_time = mean(time_to_build_model + time_to_solve_model),
  sd_time = sd(time_to_build_model + time_to_solve_model)
)
exp_7i <- full_join(t6_summ, t8_7i_summ, by = "instfname")
exp_7i_xtable <- xtable(exp_7i)
digits(exp_7i_xtable) <- c(0, 0, 0, 0, 2, 2, 0, 0, 2, 2)
print(exp_7i_xtable, only.contents = T, include.rownames = F, include.colnames = F)

# * big table at appendix with values for every instance:
#   + #var, #plates, -x% vars by round2disc, -x% plates by round2disc, (build & solve)|total time, finished, best value, best bound

t8 <- t %>% filter(experiment == 8)
t8_summ <- t8 %>% group_by(instfname) %>% summarize(
  num_vars = dplyr::first(num_vars),
  num_plates = dplyr::first(num_plates),
  mean_build_time = mean(time_to_build_model),
  mean_solve_time = mean(time_to_solve_model, na.rm = T),
  sd_solve_time = sd(time_to_solve_model, na.rm = T),
  #avg_lb = mean(obj_value),
  max_lb = max(obj_value),
  #avg_ub = mean(obj_bound),
  max_ub = min(obj_bound),
  finished = sum(stop_code == 1)
)
t8_xtable <- xtable(t8_summ)
#digits(t8_xtable) <- c(0, 0, 0, 0, 2, 2, 2, 0, 2, 0, 0)
digits(t8_xtable) <- c(0, 0, 0, 0, 2, 2, 2, 0, 0, 0)
print(t8_xtable, only.contents = T, include.rownames = F, include.colnames = F)

#sdes <- sde %>% 
#  group_by(experiment, instfname) %>% summarise(
#    mean_time_to_build_model = mean(time_to_build_model, na.rm = T),
#    sd_time_to_build_model = sd(time_to_build_model, na.rm = T),
#    mean_time_to_solve_model = mean(time_to_solve_model, na.rm = T),
#    sd_time_to_solve_model = sd(time_to_solve_model, na.rm = T),
#    finished = sum(stop_code == 1)
#  )

## transform percentage in absolute value
#tdf13 <- filter(t, experiment == 13) %>% arrange(instfname)
#tdf20 <- filter(t, experiment == 20) %>% arrange(instfname)
#tdf20$num_vars <- tdf13$num_vars * (tdf20$num_vars/100.0)
#t2 <- tdf20 %>% select(instfname, experiment, num_vars)
#t2$num_vars <- round(t2$num_vars)
#write.table(t2, file = "csv/dimitri_furini_priced_abs.csv", quote = F, sep = ";")
