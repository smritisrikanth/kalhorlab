plot_tb = tibble(Condition = c('0.4 unadjusted nonzero average',
                               '0.4 unadjusted total average',
                               '0.4 adjusted nonzero average',
                               '0.4 adjusted total average',
                               '0.65 unadjusted nonzero average',
                               '0.65 unadjusted total average',
                               '0.65 adjusted nonzero average',
                               '0.65 adjusted total average'),
                 Error = c(0.3067725, 0.2641398, 0.1364936, 0.09160725,
                           0.4363294, 0.430227, 0.1717421, 0.1477291))
plot_tb = plot_tb[order(plot_tb$Error),]
ggplot(data = plot_tb, aes(x = reorder(Condition, -Error), y = Error)) + geom_bar(stat = 'identity') +
  theme(axis.text=element_text(size=6))
