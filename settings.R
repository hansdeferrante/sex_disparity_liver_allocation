# Settings for simulations
NUMBER_YEARS_SIMULATED <- 10
TIME_HORIZON <- NUMBER_YEARS_SIMULATED*365
NUMBER_OF_CANDIDATES_PER_YEAR <- 1700
TOTAL_NUMBER_OF_CANDIDATES <- NUMBER_OF_CANDIDATES_PER_YEAR*NUMBER_YEARS_SIMULATED
INITIAL_WAITLIST_SIZE = .60*NUMBER_OF_CANDIDATES_PER_YEAR

# Distributions of donor / candidate country over the countries
D_COUNTRY_RATES <- c(G = 0.455, N = 0.135, B = 0.141, A = 0.105, H = 0.087, C = 0.056, S = 0.017)
P_COUNTRY_RATES <- c(G = 0.659, N = 0.070, B = 0.095, A = 0.040, H = 0.066, C = 0.056, S = 0.013)

# Blood type frequencies and compatibilities
BLOOD_TYPE_FREQUENCIES <- c(A = .43, O = .40, B = .12, AB = .05)
BLOOD_TYPE_COMPATIBILITIES <- list(
  A = c('A', 'AB'),
  B = c('B', 'AB'),
  AB = c('AB'),
  O = 'O'
)

# Empirical distribution of MELD at listing in ET.
EMPIRICAL_DISTRIBUTION_MELD <- c(
  `6` = 0.047, `7` = 0.035, `8` = 0.081, `9` = 0.061, `10` = 0.061, `11` = 0.06,
  `12` = 0.062, `13` = 0.064, `14` = 0.063, `15` = 0.056, `16` = 0.055, 
  `17` = 0.049, `18` = 0.042, `19` = 0.034, `20` = 0.029, `21` = 0.027, 
  `22` = 0.022, `23` = 0.021, `24` = 0.017, `25` = 0.012, `26` = 0.011, 
  `27` = 0.008, `28` = 0.008, `29` = 0.007, `30` = 0.006, `31` = 0.007, 
  `32` = 0.006, `33` = 0.005, `34` = 0.005, `35` = 0.005, `36` = 0.005, 
  `37` = 0.005, `38` = 0.004, `39` = 0.003, `40` = 0.017
)



# Set the baseline hazard function, simulate MELD scores according to the
# empirical distribution of MELD in ET. Note that MELD is assumed to be
# independent of sex. Then simulate data based on a Cox proportional hazards
# model. Parameters for this Cox proportional hazards model were obtained
# by estimating such a model on ET registry data.
f_over_c <- function(c) {function(x) pmax(0, x-c)}
f_under_c <- function(c) {function(x) pmax(0, c-x)}
f_over_20 <- f_over_c(20)
f_over_30 <- f_over_c(30)
