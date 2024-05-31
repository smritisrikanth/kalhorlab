library(ComplexHeatmap)
library(qfm)
library(tidyverse)
library(furrr)

setwd('~/Documents/QFM')
tree_panel = readRDS("../QFM/data/example/tree_panel.rds")
type_graph = tree_panel$type_graph[[5]]

mut_p = readRDS("../QFM/data/example/mut_p_marc1.rds")
#replicating each locus
mut_p$mut_rate = rep(mut_p$mut_rate, times = rep(2,50))
mut_p$recur_prob = rep(mut_p$recur_prob, times=2)
mut_p$recur_vec_list = rep(mut_p$recur_vec_list, times = rep(2,50))

#simulate data
sim_data = simulate_sc_data_mod2(type_graph, mut_p = mut_p, sample_size = 100)
saveRDS(sim_data, file = "./intermediate_data/qfm1_sim_data.rds")
sim_data = readRDS("./intermediate_data/qfm1_sim_data.rds")
#reread mutation parameters and replicate everything except mutation rate, which is the same between sites

qfm::plot_barcodes(sim_data$sc, sim_data$tr)

#get matrix of pairwise distances
dist_mat = cophenetic(sim_data$tr)


#get barcodes from two arbitrary cells
b1 = list(sim_data$sc[1,seq(1,99,2)],
          sim_data$sc[1,seq(2,100,2)])
b2 = list(sim_data$sc[24,seq(1,99,2)],
          sim_data$sc[24,seq(2,100,2)])
total_coal_time = 10.9

#diploid estimation function
estimate_coal_time <- function(b1, b2, mut_p, total_coal_time) {

        #must be adjusted to look at each item within b1 and b2
        avail_indices = !is.na(b1) & !is.na(b2)
        if (!any(avail_indices)) return(total_coal_time)
        b1 = b1[avail_indices]
        b2 = b2[avail_indices]
        mut_p = list(mut_rate = mut_p$mut_rate[avail_indices],
                     recur_vec_list = mut_p$recur_vec_list[avail_indices])

        #a = cell1 1 site 1, b = cell 1 site 2, c = cell 2 site 1, d = cell 2 site 2
        a = b1[[1]]
        b = b1[[2]]
        c = b2[[1]]
        d = b2[[2]]

        nz_ind = a != "0" & b!= "0" & c != "0" & d != "0"

        #case 1 = AB AB, case 2 = AB AC, case 3 = AB CD
        case_ind1 = (a == c & b == d) | (a == d & b == c)
        case_ind2 = (a == c | a == d | b == d | b == c) & !case_ind1
        case_ind3 = !(case_ind1 | case_ind2)


        ((c == a) | (c == d) | (d == b) | (d == a))[case_ind3]

        #case2a = AB AC
        # case 2b = AB CA
        # 2c AB CB
        #2d AB BC
        case2a = case_ind2 & (a == c)
        case2b = case_ind2 & (a == d) & !case2a
        case2c = case_ind2 & (b == d) & !case2a & !case2b
        case2d = case_ind2 & (b == c) & !case2a & !case2b & !case2c

        # case_ind1 = case_ind1 & nz_ind
        # case2a = case2a & nz_ind
        # case2b = case2b & nz_ind
        # case2c = case2c & nz_ind
        # case2d = case2d & nz_ind
        # case_ind3 = case_ind3 & nz_ind

        # assertthat::assert_that(all(case_ind1 + case2a + case2b + case2c + case2d + case_ind3 == 1))
        a1_vec = map2_dbl(mut_p$recur_vec_list[seq(1,99,2)][case_ind1], b1[[1]][case_ind1], function(r_vec, m) {
                r_vec[m]
        })
        a2a_vec = map2_dbl(mut_p$recur_vec_list[seq(1,99,2)][case2a], b1[[1]][case2a], function(r_vec, m) {
                r_vec[m]
        })
        a2b_vec = map2_dbl(mut_p$recur_vec_list[seq(1,99,2)][case2b], b1[[1]][case2b], function(r_vec, m) {
                r_vec[m]
        })
        a2c_vec = map2_dbl(mut_p$recur_vec_list[seq(1,99,2)][case2c], b1[[1]][case2c], function(r_vec, m) {
                r_vec[m]
        })
        a2d_vec = map2_dbl(mut_p$recur_vec_list[seq(1,99,2)][case2d], b1[[1]][case2d], function(r_vec, m) {
                r_vec[m]
        })
        a3_vec = map2_dbl(mut_p$recur_vec_list[seq(1,99,2)][case_ind3], b1[[1]][case_ind3], function(r_vec, m) {
                r_vec[m]
        })

        b1_vec = map2_dbl(mut_p$recur_vec_list[seq(2,100,2)][case_ind1], b1[[2]][case_ind1], function(r_vec, m) {
                r_vec[m]
        })
        b2a_vec = map2_dbl(mut_p$recur_vec_list[seq(2,100,2)][case2a], b1[[2]][case2a], function(r_vec, m) {
                r_vec[m]
        })
        b2b_vec = map2_dbl(mut_p$recur_vec_list[seq(2,100,2)][case2b], b1[[2]][case2b], function(r_vec, m) {
                r_vec[m]
        })
        b2c_vec = map2_dbl(mut_p$recur_vec_list[seq(2,100,2)][case2c], b1[[2]][case2c], function(r_vec, m) {
                r_vec[m]
        })
        b2d_vec = map2_dbl(mut_p$recur_vec_list[seq(2,100,2)][case2d], b1[[2]][case2d], function(r_vec, m) {
                r_vec[m]
        })
        b3_vec = map2_dbl(mut_p$recur_vec_list[seq(2,100,2)][case_ind3], b1[[2]][case_ind3], function(r_vec, m) {
                r_vec[m]
        })


        c3_vec = map2_dbl(mut_p$recur_vec_list[seq(1,99,2)][case_ind3], b2[[1]][case_ind3], function(r_vec, m) {
                r_vec[m]
        })
        d_vec = map2_dbl(mut_p$recur_vec_list[seq(2,100,2)][case_ind3], b2[[2]][case_ind3], function(r_vec, m) {
                r_vec[m]
        })

        c2a_vec = map2_dbl(mut_p$recur_vec_list[seq(2,100,2)][case2a], b2[[2]][case2a], function(r_vec, m) {
                r_vec[m]
        })
        c2b_vec = map2_dbl(mut_p$recur_vec_list[seq(1,99,2)][case2b], b2[[1]][case2b], function(r_vec, m) {
                r_vec[m]
        })
        c2c_vec = map2_dbl(mut_p$recur_vec_list[seq(1,99,2)][case2c], b2[[1]][case2c], function(r_vec, m) {
                r_vec[m]
        })
        c2d_vec = map2_dbl(mut_p$recur_vec_list[seq(2,100,2)][case2d], b2[[2]][case2d], function(r_vec, m) {
                r_vec[m]
        })

        est_func <- function(t_val) {
                t_val = seq(0,11,0.1)
                y_val <- map_dbl(seq(0,11,0.1), function(t_val) {
                        e1 = est_func_case3(t = t_val,
                                            lambda = mut_p$mut_rate[case_ind1],
                                            t_total = 10.9,
                                            a = a1_vec,
                                            b = b1_vec
                        )
                        e2a = est_func_case4(t = t_val,
                                             lambda = mut_p$mut_rate[case2a],
                                             t_total = 10.9,
                                             a = a2a_vec,
                                             b = b2a_vec,
                                             c = c2a_vec
                        )
                        e2b = est_func_case4(t = t_val,
                                             lambda = mut_p$mut_rate[case2b],
                                             t_total = 10.9,
                                             a = a2b_vec,
                                             b = b2b_vec,
                                             c = c2b_vec
                        )
                        e2c = est_func_case4(t = t_val,
                                             lambda = mut_p$mut_rate[case2c],
                                             t_total = 10.9,
                                             a = b2c_vec,
                                             b = a2c_vec,
                                             c = c2c_vec
                        )
                        e2d = est_func_case4(t = t_val,
                                             lambda = mut_p$mut_rate[case2d],
                                             t_total = 10.9,
                                             a = b2d_vec,
                                             b = a2d_vec,
                                             c = c2d_vec
                        )
                        e3 = est_func_case5(t = t_val,
                                            lambda = mut_p$mut_rate[case_ind3],
                                            t_total = 10.9,
                                            a = a3_vec,
                                            b = b3_vec,
                                            c = c3_vec,
                                            d = d_vec
                        )
                        # sum(c(e1, e2a, e2b, e2c, e2d, e3))
                        sum(c(e3))
                })


                plot(t_val,y_val,type = 'l')
                plot(t,c(e1, e2a, e2b, e2c, e2d, e3))
                abline(h=0)
                # sum(c(e1, e2a, e2b, e2c, e2d, e3))
                sum(c(e1, e2a, e2b, e2c, e2d, e3))
        }
        if (est_func(10.9) > 0 & est_func(0.1) > 0) {
                return(total_coal_time)
        }
        total_coal_time = 10.9
        if (est_func(10.8) < 0 & est_func(0.1) < 0) {
                return(0)
        }
        t_est = uniroot(est_func, interval = c(0.1, 10.9))$root
        t_est
}

#returns d(log(p))/dt for case 1
est_func_case3 <- function(t, lambda, t_total, a, b) {
        e_lt = exp(-1 * lambda * t)
        e_lt_total = exp(-1 * lambda * t_total)
        e_lt_total_minus = 1 - exp(-1 * lambda * (t_total-t))

        a_na = is.na(a)
        b_na = is.na(b)

        lambda_a = lambda[a_na]
        lambda_b = lambda[b_na]

        e_lt_total_minus_a = 1 - exp(-1 * lambda_a * (t_total-t))
        e_lt_total_minus_b = 1 - exp(-1 * lambda_b * (t_total-t))

        #inserts correct probability in every spot in the allele vector that was originally NA (meaning the corresponding allele was '0')
        a[is.na(a)] = (e_lt_total_minus_a + 1)/e_lt_total_minus_a
        a[is.infinite(a)] = 1

        b[is.na(b)] = (e_lt_total_minus_b + 1)/e_lt_total_minus_b
        b[is.infinite(b)] = 1

        term_1 = 2*a*b*lambda*e_lt*(1-e_lt) +
                (b*a^2 + a*b^2)*(-1*lambda*e_lt*e_lt_total_minus^3 - 3*lambda*e_lt_total*e_lt_total_minus^2) +
                (a^2 * b^2) * (-6*lambda*e_lt^2*e_lt_total_minus^4 - 12*lambda*(e_lt_total_minus + 1)*e_lt_total_minus^3)
        term_2 = a * b * (1 - e_lt)^2 + e_lt * (a^2*b + a*b^2) * (e_lt_total_minus)^3 + 3 * e_lt^2 * a^2 * b^2 * (e_lt_total_minus)^4
        term_2[term_2 == 0] = 1
        term_1[term_2 == 0] = 0
        term_1/term_2
}

#returns d(log(p))/dt for case 2
est_func_case4 <- function(t, lambda, t_total, a, b, c) {

        e_lt = exp(-1 * lambda * t)
        e_lt_total = exp(-1 * lambda * t_total)
        e_lt_total_minus = 1 - exp(-1 * lambda * (t_total-t))
        e_lt_total_plus = exp(-1 * lambda * (t_total+t))

        a_na = is.na(a)
        b_na = is.na(b)
        c_na = is.na(c)

        lambda_a = lambda[a_na]
        lambda_b = lambda[b_na]
        lambda_c = lambda[c_na]

        e_lt_total_minus_a = 1 - exp(-1 * lambda_a * (t_total-t))
        e_lt_total_minus_b = 1 - exp(-1 * lambda_b * (t_total-t))
        e_lt_total_minus_c = 1 - exp(-1 * lambda_c * (t_total-t))

        a[is.na(a)] = (e_lt_total_minus_a + 1)/e_lt_total_minus_a
        b[is.na(b)] = (e_lt_total_minus_b + 1)/e_lt_total_minus_b
        c[is.na(c)] = (e_lt_total_minus_c + 1)/e_lt_total_minus_c

        a[is.infinite(a)] = 1
        b[is.infinite(b)] = 1
        c[is.infinite(c)] = 1

        term_1 = a*b*c*((-2*lambda*e_lt*e_lt_total_minus^2 - 4*lambda*e_lt_total*e_lt_total_minus)*(1 - e_lt + a*e_lt*e_lt_total_minus^2)) +
                a*b*c*-2*lambda*e_lt*e_lt_total_minus^2*(lambda*e_lt - a*lambda*e_lt*e_lt_total_minus^2 - 2*a*e_lt_total*e_lt_total_minus) +
                a^2*b*c*(-8*lambda*e_lt^2*e_lt_total_minus^4 - 16*lambda*e_lt_total_plus*e_lt_total_minus^3)
        term_2 = 2 * e_lt * b * c * e_lt_total_minus^2 * (a * (1 - e_lt) + a^2 * e_lt * e_lt_total_minus^2) +
                4 * a^2 * b * c * e_lt^2 * e_lt_total_minus^4
        term_1[term_2 == 0] = 0
        term_2[term_2 == 0] = 1
        term_1/term_2
}

est_func_case5 <- function(t, lambda, t_total, a, b, c, d) {

        # a = a3_vec
        # b = b3_vec
        # c = c3_vec
        # d = d_vec
        #
        # t = 1.7
        # t_total = 10.9

        e_lt_total_minus = exp(-1 * lambda * (t_total-t))
        e_lt = exp(-1 * lambda * t)

        a_na = is.na(a)
        b_na = is.na(b)
        c_na = is.na(c)
        d_na = is.na(d)

        lambda_a = lambda[a_na]
        lambda_b = lambda[b_na]
        lambda_c = lambda[c_na]
        lambda_d = lambda[d_na]

        e_lt_total_minus_a = exp(-1 * lambda_a * (t_total-t))
        e_lt_total_minus_b = exp(-1 * lambda_b * (t_total-t))
        e_lt_total_minus_c = exp(-1 * lambda_c * (t_total-t))
        e_lt_total_minus_d = exp(-1 * lambda_d * (t_total-t))

        a[is.na(a)] = (e_lt_total_minus_a)/(1-e_lt_total_minus_a)
        b[is.na(b)] = (e_lt_total_minus_b)/(1-e_lt_total_minus_b)
        c[is.na(c)] = (e_lt_total_minus_c)/(1-e_lt_total_minus_c)
        d[is.na(d)] = (e_lt_total_minus_d)/(1-e_lt_total_minus_d)

        a[is.infinite(a)] = 1
        b[is.infinite(b)] = 1
        c[is.infinite(c)] = 1
        d[is.infinite(d)] = 1

        term_1 = -2*lambda - (4*lambda*e_lt_total_minus)/(1 - e_lt_total_minus)
        #term_2 = 4 * e_lt^2 * a * b * c * d * e_lt_total_minus^4

        term_1[is.infinite(term_1)] = -2*lambda
        #term_1[term_2 == 0] = 0
        #term_2[term_2 == 0] = 1
        term_1
}

#testing
b1 = list(sim_data$sc[100,seq(1,99,2)],
          sim_data$sc[100,seq(2,100,2)])
b2 = list(sim_data$sc[10,seq(1,99,2)],
          sim_data$sc[10,seq(2,100,2)])

est_t <- estimate_coal_time(b1, b2, mut_p, 10.9)
dist_mat["type_-1_gen_4_1", "type_-2_gen_4_1"]


#original phylotime
qfm::est_func_case1 <- function (t, lambda, a, t_total) 
{
  e_lt = exp(lambda * t)
  e_lt_total = exp(lambda * t_total)
  lambda * ((a - 1) * e_lt^2 - a)/((e_lt_total/e_lt + a - 1) * 
                                     e_lt^2 - 2 * a * e_lt + a)
}

qfm::est_func_case1a <- function (t, lambda, t_total) 
{
  e_lt = exp(lambda * t)
  e_lt_total = exp(lambda * t_total)
  lambda/(1 - e_lt_total/e_lt)
}

qfm::est_func_case2 <- function (t, i_total, lambda) 
{
  e_lt_inv = exp(lambda * t)
  lambda * ((i_total - 1) * e_lt_inv + 1)/(e_lt_inv - 1)
}

estimate_coal_time <- function (b1, b2, mut_p, total_coal_time, min_alele_prob = 0) 
{
  avail_indices = !is.na(b1) & !is.na(b2)
  if (!any(avail_indices)) 
    return(total_coal_time)
  b1 = b1[avail_indices]
  b2 = b2[avail_indices]
  mut_p = list(mut_rate = mut_p$mut_rate[avail_indices], recur_vec_list = mut_p$recur_vec_list[avail_indices])
  case_ind = (b1 == b2 & b1 != "0" & b2 != "0")
  a_vec = map2_dbl(mut_p$recur_vec_list[case_ind], b1[case_ind], 
                   function(r_vec, m) {
                     r_vec[m]
                   })
  a_ind = a_vec > min_alele_prob
  est_func <- function(t_val) {
    c1 = est_func_case1(t = t_val, lambda = mut_p$mut_rate[case_ind][a_ind], 
                        a = a_vec[a_ind], t_total = total_coal_time)
    c1a = est_func_case1a(t = t_val, lambda = mut_p$mut_rate[case_ind][!a_ind], 
                          t_total = total_coal_time)
    c2 = est_func_case2(t = t_val, i_total = (b1[!case_ind] != 
                                                "0") + (b2[!case_ind] != "0"), lambda = mut_p$mut_rate[!case_ind])
    sum(c(c1, c1a, c2))
  }
  if (est_func(total_coal_time) > 0 & est_func(0.1) > 0) {
    return(total_coal_time)
  }
  if (est_func(total_coal_time) < 0 & est_func(0.1) < 0) {
    return(0)
  }
  t_est = uniroot(est_func, interval = c(0.1, total_coal_time))$root
  t_est
}

