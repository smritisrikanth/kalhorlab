num_steps = 800
start_val = 1
delta = 10^-10
eta_vals = c(0.2,0.8)
v0_vals = c(0.01,0.1)
esd_vals = c(step_size/10, step_size/20, step_size/100)
num_rep = 5

plot_tb = data.frame()
for (eta in eta_vals) {
  for (v0 in v0_vals) {
    for (esd in esd_vals) {
      for (j in 1:num_rep) {
        error = rnorm(num_steps-1,0,esd)
        coefs = cumsum(map_dbl(0:(num_steps-1), function(i) eta^i))
        error_terms = c(0,map_dbl(1:length(error), function(i) sum(error[1:i]*rev(coefs[1:i]))))
        walk = c(start_val,coefs*v0 + error_terms + start_val)
        w_stable = start_val + v0/(1-eta)
        t_stable = -log((1-eta)/delta+1)/log(eta) * step_size
        plot_tb = rbind(plot_tb,
                        tibble(time = linspace(0,8,num_steps+1),
                               value = walk,
                               eta = eta,
                               v0 = v0,
                               esd = esd,
                               w_stable = w_stable,
                               t_stable = t_stable,
                               rep = j))
        
      }
    }
  }
}

plot = ggplot(data = plot_tb, aes(x = time, y = value, group = rep)) +
  geom_line() + facet_wrap(~esd+v0+eta, ncol = 4) +
  geom_vline(aes(xintercept = t_stable),
             linetype = 2,
             color = 'red') +
  geom_hline(aes(yintercept = w_stable),
           linetype = 2,
           color = 'blue') +
  xlab('Time (Days)') + ylab('Value (AU)') +
  theme_pubr() + theme(legend.position = 'right')

png(filename = './eta_v0_figure.png',
    width = 5, height = 5, units = "in", res = 1000)
print(plot)
dev.off()

# eta = 0.1
# v0 = 0.1
# 
# error = rnorm(num_steps-1,0,step_size)
# coefs = cumsum(map_dbl(0:(num_steps-1), function(i) eta^i))
# error_terms = c(0,error * coefs[1:length(error)])
# walk = c(start_val,coefs*v0 + error_terms + start_val)
# 
# plot_tb = tibble(time = linspace(0,num_steps*step_size,num_steps+1),
#                  value = walk)
# plot_tb$eta = eta
# plot_tb1 = plot_tb
# 
# plot1 = ggplot(data = plot_tb, aes(x = time, y = value)) +
#   geom_line() +
#   xlab('Time (Days)') + ylab('Value (AU)') +
#   ggtitle(paste0('Eta = ', eta, ', v0 = ', v0))
# plot1
# 
# eta = 0.5
# v0 = 0.1
# 
# error = rnorm(num_steps-1,0,step_size)
# coefs = cumsum(map_dbl(0:(num_steps-1), function(i) eta^i))
# error_terms = c(0,error * coefs[1:length(error)])
# walk = c(start_val,coefs*v0 + error_terms + start_val)
# 
# plot_tb = tibble(time = linspace(0,num_steps*step_size,num_steps+1),
#                  value = walk)
# plot_tb$eta = eta
# plot_tb2 = plot_tb
# 
# plot2 = ggplot(data = plot_tb, aes(x = time, y = value)) +
#   geom_line() +
#   xlab('Time (Days)') + ylab('Value (AU)') +
#   ggtitle(paste0('Eta = ', eta, ', v0 = ', v0))
# plot2
# 
# eta = 0.8
# v0 = 0.1
# 
# error = rnorm(num_steps-1,0,step_size)
# coefs = cumsum(map_dbl(0:(num_steps-1), function(i) eta^i))
# error_terms = c(0,error * coefs[1:length(error)])
# walk = c(start_val,coefs*v0 + error_terms + start_val)
# 
# plot_tb = tibble(time = linspace(0,num_steps*step_size,num_steps+1),
#                  value = walk)
# plot_tb$eta = eta
# plot_tb3 = plot_tb
# 
# plot3 = ggplot(data = plot_tb, aes(x = time, y = value)) +
#   geom_line() + ylim(0,0.8) +
#   xlab('Time (Days)') + ylab('Value (AU)') +
#   ggtitle(paste0('Eta = ', eta, ', v0 = ', v0))
# plot3
# 
# eta = 0.11
# v0 = 0.01
# 
# error = rnorm(num_steps-1,0,step_size)
# coefs = cumsum(map_dbl(0:(num_steps-1), function(i) eta^i))
# error_terms = c(0,error * coefs[1:length(error)])
# walk = c(start_val,coefs*v0 + error_terms + start_val)
# 
# plot_tb = tibble(time = linspace(0,num_steps*step_size,num_steps+1),
#                  value = walk)
# plot_tb$eta = eta
# plot_tb4 = plot_tb
# 
# plot4 = ggplot(data = plot_tb, aes(x = time, y = value)) +
#   geom_line() + ylim(0,0.8) +
#   xlab('Time (Days)') + ylab('Value (AU)') +
#   ggtitle(paste0('Eta = ', eta, ', v0 = ', v0))
# plot4
# 
# 
# (plot1 + plot2) / (plot3 + plot4)
# 
# 
# plot_tb = rbind(plot_tb1,plot_tb2, plot_tb3, plot_tb4)
# plot = ggplot(data = plot_tb, aes(x = time, y = value)) +
#   geom_line() + facet_wrap(~eta) +
#   xlab('Time (Days)') + ylab('Value (AU)') +
#   theme_pubr()
# 
# #random walk plots
# eta = 0.8
# v0 = 0.1
# start_val = 1
# num_steps = 6/step_size
# 
# pdf(file = paste0('eta0.8_v0_0.1_errorstd0.001.pdf'))
# plot_tb = data.frame()
# for (i in 1:10) {
#   error = rnorm(num_steps-1,0,step_size/100)
#   coefs = cumsum(map_dbl(0:(num_steps-1), function(i) eta^i))
#   error_terms = c(0,map_dbl(1:length(error), function(i) sum(error[1:i]*rev(coefs[1:i]))))
#   walk = c(start_val,coefs*v0 + error_terms + start_val)
#   
#   # walk = rep(0,num_steps+1)
#   # walk[1] = start_val
#   # velocity = rep(0,num_steps+1)
#   # velocity[1] = v0
#   # for (t in 1:(num_steps)) {
#   #   walk[t+1] = walk[t] + velocity[t]
#   #   e = rnorm(1,0,step_size)
#   #   velocity[t+1] = eta*(velocity[t] + e)
#   # }
#   
#   x = tibble(time = linspace(0,6,num_steps+1),
#              value = walk,
#              rep = rep(i,num_steps+1))
#   plot_tb = rbind(plot_tb, x)
# }
# dev.off()
# 
# plot_tb$rep = factor(plot_tb$rep)
# 
# ggplot(data = plot_tb, aes(x = time, y = value, color = rep)) +
#   geom_line()
# 
