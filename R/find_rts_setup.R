get_ps <- function(si_shape, si_rate, k = 100){
  # might want to rename this to get_si (since that's what it's doing...)
  #
  #' Find Discrete Serial Interval
  #'
  #' Gets discrete probability distribution for serial interval using given
  #'   shape and rate parameters.
  #'
  #' @param si_shape A float, shape parameter for gamma distribution.
  #' @param si_rate A float, rate parameter for gamma distribution.
  #' @param k An integer, maximum length of gamma distribution
  #' @return List with values for discrete gamma distribution with length
  #'   \code{k}.
  #' @examples getPs(si_shape=1.47, si_rate=0.04529, k=100)

  ps_si <- rep(0, k)
  for (i in 1:k){
    ps_si[i] <- pgamma(i, shape = si_shape, rate = si_rate) -
                pgamma(i - 1, shape = si_shape, rate = si_rate)
  }

  # make probabilities sum to 1
  ps_si <- ps_si / sum(ps_si)

  return(ps_si)
}


get_si_weights <- function(t_list, ps_si){
  # can change this to take in time column, convert to weeknum, etc
  #
  #' Serial Interval Based Probability of Infection
  #'
  #' Gets probability of infection based solely on serial interval
  #'   (ie. time between primary and secondary infection.)
  #'
  #' @param t_list List of times. For weekly data, this list will be the week
  #'   number starting at 1.
  #' @param ps_si List of discrete serial interval probabilities. Result of
  #'   \code{get_ps}.
  #' @return Matrix of probabilities of infection based solely on serial
  #'   interval.
  #' @example get_si_weights(seq(1:200), get_ps(si_shape=1.47, si_rate=0.04529,
  #'   k=100))

  si_weights <- matrix(nrow = length(t_list), ncol = length(t_list))
  for (t1 in 2:length(t_list)){
    for (t2 in 1:length(t_list)){
      date_diff <- t_list[t1] - t_list[t2]
      p_si <- ifelse(date_diff > 0, ps_si[date_diff], 0)
      si_weights[t1, t2] <- p_si
    }
  }
  si_weights[1, ] <- 0
  si_weights[is.na(si_weights)] <- 0

  return(si_weights)
}


get_rts <- function(data, si_weights){
  #' Finds Reproductive Numbers
  #'
  #' Finds reproductive numbers based on number of individuals infected at each
  #'   time point. Reweighs the \code{si_weights} by number infected.
  #'
  #' @param data Dataframe to be used for predicting reproductive number.
  #' @param si_weights Matrix of probabilities of infection based on SI. Output
  #'   from \code{get_si_weights}.
  #' @return List of effective reproductive numbers based on the time list
  #'   passed into \code{get_si_weights}.
  #' @section Warning:
  #' Data must contain 'cases' column.

  reweighed_row_sum <- si_weights %*% data$cases
  reweighed_prob <- si_weights / as.vector(reweighed_row_sum)
  reweighed_prob[is.na(reweighed_prob)] <- 0

  rts <- colSums(reweighed_prob * data$cases)

  return(rts)
}
