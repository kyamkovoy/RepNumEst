find_rts <- function(data, si_shape=1.47, si_rate=0.04529, k=100){
  #' Finds Reproductive Numbers.
  #'
  #' Using data and specified serial interval, finds the effective reproductive
  #'   number. The default serial interval uses an estimate for the United
  #'   States on the weekly time scale.
  #'
  #' @param data Dataframe to be used for estimating reproductive numbers. Must
  #'   contain 'cases' column.
  #' @param si_shape Shape parameter for serial interval.
  #' @param si_rate Rate parameter for serial interval.
  #' @param k Maximum length of serial interval.
  #' @return List of effective reproductive numbers and corresponding week
  #'   numbers.

  dates <- seq(1, nrow(data))

  # get discrete gamma distribution for serial interval
  ps_si <- rep(0, k)
  for (i in 1:k){
    ps_si[i] <- pgamma(i, shape = si_shape, rate = si_rate) -
                pgamma(i - 1, shape = si_shape, rate = si_rate)
  }
  # ensure that ps sum to 1
  ps_si <- ps_si / sum(ps_si)


  # get weights from serial distribution
  si_weights <- matrix(nrow = length(dates), ncol = length(dates))
  for (t1 in 2:length(dates)){
    for (t2 in 1:length(dates)){
      date_diff <- dates[t1] - dates[t2]
      p_si <- ifelse(date_diff > 0, ps_si[date_diff], 0)
      si_weights[t1, t2] <- p_si
    }
  }
  si_weights[1, ] <- 0
  si_weights[is.na(si_weights)] <- 0


  # get reproductive numbers
  reweighed_row_sum <- si_weights %*% data$cases
  reweighed_prob <- si_weights / as.vector(reweighed_row_sum)
  reweighed_prob[is.na(reweighed_prob)] <- 0

  rts <- colSums(reweighed_prob * data$cases)

  return(list(rts = rts, date = dates))  # might wanna rename these
  # SAVE OUTPUT
}


find_rts2 <- function(cases, si_shape=1.47, si_rate=0.04529, k=100){
  #' Finds Reproductive Numbers.
  #'
  #' Using data and specified serial interval, finds the effective reproductive
  #'   number. The default serial interval uses an estimate for the United
  #'   States on the weekly time scale. Only requires a list of cases and not
  #'   the entire dataframe for the data.
  #'
  #' @param cases List of cases to be used for estimating reproductive numbers.
  #' @param si_shape Shape parameter for serial interval.
  #' @param si_rate Rate parameter for serial interval.
  #' @param k Maximum length of serial interval.
  #' @return List of effective reproductive numbers and corresponding week
  #'   numbers.

  dates <- seq(1, length(cases))

  # get discrete gamma distribution density
  ps_si <- rep(0, k)
  for (i in 1:k){
    ps_si[i] <- pgamma(i, shape = si_shape, rate = si_rate) -
                pgamma(i - 1, shape = si_shape, rate = si_rate)
  }
  # ensure probabilities sum to 1
  ps_si <- ps_si / sum(ps_si)

  # get weights from serial distribution
  si_weights <- matrix(nrow = length(dates), ncol = length(dates))
  for (t1 in 2:length(dates)){
    for (t2 in 1:length(dates)){
      date_diff <- dates[t1] - dates[t2]
      p_si <- ifelse(date_diff > 0, ps_si[date_diff], 0)
      si_weights[t1, t2] <- p_si
    }
  }
  si_weights[1, ] <- 0
  si_weights[is.na(si_weights)] <- 0


  # get reproductive numbers
  reweighed_row_sum <- si_weights %*% cases
  reweighed_prob <- si_weights / as.vector(reweighed_row_sum)
  reweighed_prob[is.na(reweighed_prob)] <- 0

  rts <- colSums(reweighed_prob * cases)

  return(list(rts = rts, date = dates))  # might want to rename these
  # SAVE OUTPUT
}

find_rts_region <- function(data, reg_mat, si_shape=1.47, si_rate=0.04529,
                         k=100){
  #' Finds Region Specific Reproductive Numbers
  #'
  #' Using data and specified serial interval, finds the effective reproductive
  #'   for each region while allowing for interaction between regions. The
  #'   default serial interval uses an estimate for the United States on the
  #'   weekly time scale.
  #'
  #' @param data Dataframe to be used for estimating reproductive numbers. Must
  #'   contain 'cases' column and numerical 'region' column.
  #' @param reg_mat Matrix (4 x 4) representing interactions between regions.
  #' @param si_shape Shape parameter for serial interval.
  #' @param si_rate Rate parameter for serial interval.
  #' @param k Maximum length of serial interval.
  #' @return List of effective reproductive numbers for each regions and
  #'   corresponding week numbers.

  dates <- seq(1, length(factor(data$epiweek)) / 4)

  # get discrete gamma distribution density
  ps_si <- rep(0, k)
  for (i in 1:k){
    ps_si[i] <- pgamma(i, shape = si_shape, rate = si_rate) -
                pgamma(i - 1, shape = si_shape, rate = si_rate)
  }

  ps_si <- ps_si / sum(ps_si)

  # get weights from serial distribution
  si_weights <- matrix(nrow = length(dates), ncol = length(dates))
  for (t1 in 2:length(dates)){
    for (t2 in 1:length(dates)){
      date_diff <- dates[t1] - dates[t2]
      p_si <- ifelse(date_diff > 0, ps_si[date_diff], 0)
      si_weights[t1, t2] <- p_si
    }
  }
  si_weights[1, ] <- 0
  si_weights[is.na(si_weights)] <- 0

  # reweigh these weights by the interaction between regions
  int_weights <- si_weights %x% reg_mat

  # list of all cases, staggered by region
  cases <- data$cases
  cases_0s <- cases # need to 0 out the entries where there were no observations
  cases_0s[cases > 0] <- 1

  int_weights2 <- t(t(int_weights) * cases_0s)

  # weigh by the number of cases
  reweighed_row_sum <- int_weights2 %*% cases
  reweighed_prob <- int_weights2 / as.vector(reweighed_row_sum)
  reweighed_prob[is.na(reweighed_prob)] <- 0

  rts_raw <- colSums(reweighed_prob * cases)  # calculate reproductive numbers

  # get reproductive numbers per region
  temp <- nrow(data)
  rts1 <- rts_raw[seq(1, temp, 4)]
  rts2 <- rts_raw[seq(2, temp, 4)]
  rts3 <- rts_raw[seq(3, temp, 4)]
  rts4 <- rts_raw[seq(4, temp, 4)]

  return(list(rts_mw = rts1, rts_ne = rts2, rts_s = rts3, rts_w = rts4,
              week = dates))
}


find_rts_strata <- function(data, num_strata=4, strata_mat=diag(4),
                          si_shape=1.47, si_rate=0.04529, k=100){

  #' Finds Region Specific Reproductive Numbers
  #'
  #' Using data and specified serial interval, finds the effective reproductive
  #'   for each region while allowing for interaction between regions. The
  #'   default serial interval uses an estimate for the United States on the
  #'   weekly time scale.
  #'
  #' @param data Dataframe to be used for estimating reproductive numbers. Must
  #'   contain 'cases' column and numerical 'region' column.
  #' @param num_strata Integer for the number of strata that the data is being
  #'   split by. Default is 4.
  #' @param strata_mat Matrix representing interactions between strata.
  #' @param si_shape Shape parameter for serial interval.
  #' @param si_rate Rate parameter for serial interval.
  #' @param k Maximum length of serial interval.
  #' @return List of effective reproductive numbers for each stratum and
  #'   corresponding week numbers. For analysis, these reproductive numbers must
  #'   be seperated.

  dates <- seq(1, nrow(data) / num_strata)

  # get discrete gamma distribution density
  ps_si <- rep(0, k)
  for (i in 1:k){
    ps_si[i] <- pgamma(i, shape = si_shape, rate = si_rate) -
                pgamma(i - 1, shape = si_shape, rate = si_rate)
  }
  ps_si <- ps_si / sum(ps_si)


  # get weights from serial distribution
  si_weights <- matrix(nrow = length(dates), ncol = length(dates))
  for (t1 in 2:length(dates)){
    for (t2 in 1:length(dates)){
      date_diff <- dates[t1] - dates[t2]
      p_si <- ifelse(date_diff > 0, ps_si[date_diff], 0)
      si_weights[t1, t2] <- p_si
    }
  }
  si_weights[1, ] <- 0
  si_weights[is.na(si_weights)] <- 0

  # reweigh these weights by the interaction between strata
  int_weights <- si_weights %x% strata_mat

  # list of all cases, staggered by stratum
  cases <- data$cases
  cases_0s <- cases # need to 0 out the entries where there were no observations
  cases_0s[cases > 0] <- 1

  int_weights2 <- t(t(int_weights) * cases_0s)

  # weigh by the number of cases
  reweighed_row_sum <- int_weights2 %*% cases
  reweighed_prob <- int_weights2 / as.vector(reweighed_row_sum)
  reweighed_prob[is.na(reweighed_prob)] <- 0

  rts_raw <- colSums(reweighed_prob * cases)  # calculate reproductive numbers

  return(rts_raw)
}
