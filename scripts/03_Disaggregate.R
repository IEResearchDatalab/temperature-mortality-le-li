################################################################################
#
# Temperature-mortality to life expectancy and lifespan inequality
#
# Pipeline: 03_Disaggregate.R
#   PCLM disaggregation of ANs from wide age groups to single years
#   Parallelised over cities x SSP x GCM
#
################################################################################

if (length(ls()) == 0) source("02_ComputeAN.R")

library(ungroup)

#---------------------------
# Load annual ANs
#---------------------------

annual_ans <- fread(file.path(path_out, "ans_annual_all_cities.csv"))

city_gcm_ssp <- unique(annual_ans[, .(city_code, gcm, ssp)])
n_combos <- nrow(city_gcm_ssp)

cat(sprintf("Disaggregating %d city-GCM-SSP combinations...\n", n_combos))

#---------------------------
# Prepare parallelisation
#---------------------------

cl <- makeCluster(ncores)
registerDoParallel(cl)

writeLines(c(""), "temp/log_pclm.txt")
cat(as.character(as.POSIXct(Sys.time())), file = "temp/log_pclm.txt", append = TRUE)

#---------------------------
# PCLM function
#---------------------------

disaggregate_an <- function(an_values, age_starts, nlast = 16,
  kr = 2, deg = 3) {

  M <- pclm(x = age_starts, y = an_values, nlast = nlast,
    out.step = 1, control = list(lambda = NA, kr = kr, deg = deg,
      max.iter = 200, tol = 1e-8, opt.method = "BIC"))
  as.numeric(fitted(M))
}

#---------------------------
# Parallel loop
#---------------------------

all_single <- foreach(i = seq_len(n_combos),
  .packages = c("data.table", "ungroup"),
  .combine = function(x, y) rbind(x, y, fill = TRUE),
  .errorhandling = "pass") %dopar% {

  cc <- city_gcm_ssp$city_code[i]
  gc <- city_gcm_ssp$gcm[i]
  ss <- city_gcm_ssp$ssp[i]

  if (i %% 50 == 0)
    cat("\n", "i = ", i, cc, gc, as.character(Sys.time()), "\n",
      file = "temp/log_pclm.txt", append = TRUE)

  city_data <- annual_ans[city_code == cc & gcm == gc & ssp == ss]

  results_list <- list()
  combos <- unique(city_data[, .(year, temp_range)])

  age_starts <- c(20, 45, 65, 75, 85)

  for (j in seq_len(nrow(combos))) {
    yr <- combos$year[j]
    tr <- combos$temp_range[j]

    subset_dt <- city_data[year == yr & temp_range == tr]
    an_vec <- subset_dt$an_est

    if (length(an_vec) != 5) next

    single_an <- tryCatch(
      disaggregate_an(an_vec, age_starts),
      error = function(e) {
        # Uniform distribution across single ages within each wide group
        rep(an_vec / c(25, 20, 10, 10, 16), times = c(25, 20, 10, 10, 16))
      })

    results_list[[length(results_list) + 1]] <- data.table(
      city_code = cc, gcm = gc, ssp = ss,
      year = yr, age = 20:100, temp_range = tr,
      AN = single_an)
  }

  if (length(results_list) > 0)
    rbindlist(results_list)
  else
    NULL
}

stopCluster(cl)

fwrite(all_single, file.path(path_out, "ans_single_age_all_cities.csv"))
cat(sprintf("Single-age AN: %s rows\n",
  format(nrow(all_single), big.mark = ",")))