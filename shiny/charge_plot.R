require(ggplot2)

charges <-read.csv("charge_curves2.csv", header = TRUE, stringsAsFactors = TRUE)

# The palette with grey:
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

#ylim(0,3) +
#  xlim(0,1) +

p <- 
  ggplot(data=charges, aes(x=CDF, y=ER, color=TYPE)) + 
  geom_line(size=1) + #aes(linetype=TYPE), size=1) +  
  scale_y_continuous(limits=c(0,3), breaks=NULL) +
  scale_x_continuous(limits=c(0,1), breaks=0:1) +
  scale_colour_manual(values=c("#999999", "#E69F00", "#56B4E9"), 
                      name="Adjustment",
                      breaks=c("M18","M30","ULT"),
                      labels=c("first","second","ultimate")) +
  xlab("cummulative probability") + ylab("aggregate loss") + 
  ggtitle("Aggregate retained layer loss over time") +
  theme_bw() +
  theme(legend.position=c(.15, .8)) +
  theme(plot.title = element_text(lineheight=.8, face="bold"))
p

p <- 
  ggplot(data=charges, aes(x=CDF, y=ER, color=TYPE)) + 
  geom_line(size=1) + #aes(linetype=TYPE), size=1) +  
  scale_y_continuous(limits=c(0,3), breaks=NULL) +
  scale_x_continuous(limits=c(0,1), breaks=0:1) +
  scale_colour_manual(values=c("#999999", "#E69F00", "#56B4E9"), 
                      name="Adjustment",
                      breaks=c("M18","M30","ULT"),
                      labels=c("first","second","ultimate")) +
  xlab("cummulative probability") + ylab("aggregate loss") + 
  ggtitle("Aggregate retained layer loss over time") +
  theme_bw() +
  theme(legend.position=c(.15, .8)) +
  theme(plot.title = element_text(lineheight=.8, face="bold"))
p


pdf("plots.pdf")
p
dev.off()



bp + scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"), 
                       name="Experimental\nCondition",
                       breaks=c("ctrl", "trt1", "trt2"),
                       labels=c("Control", "Treatment 1", "Treatment 2"))



# To use for fills, add
scale_fill_manual(values=cbPalette)

# To use for line and point colors, add
scale_colour_manual(values=cbPalette)