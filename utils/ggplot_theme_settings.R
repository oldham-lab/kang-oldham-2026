options(bitmapType = 'cairo')

library(ggplot2)

th <- theme_classic() +
  theme(plot.title=element_text(hjust=0.5, size=14, face="bold"),
        plot.subtitle=element_text(hjust=0.5,size=12),
        axis.title.x=element_text(size=12),
        axis.title.y=element_text(size=12),
        axis.text.x=element_text(size=10),
        axis.text.y=element_text(size=10),
        legend.title=element_blank(),
        legend.text=element_text(size=10),
        legend.position="bottom",
        #text=element_text(family="Arial"),
        plot.margin=margin(1,1,1,1,"cm")) 

th2 <- theme_bw() +
  theme(plot.title=element_text(hjust=0.5, size=14, face="bold"),
        plot.subtitle=element_text(hjust=0.5,size=12),
        axis.title.x=element_text(size=12),
        axis.title.y=element_text(size=12),
        axis.text.x=element_text(size=10),
        axis.text.y=element_text(size=10),
        legend.title=element_blank(),
        legend.text=element_text(size=10),
        legend.position="bottom",
        text=element_text(family="Arial"),
        plot.margin=margin(1,1,1,1,"cm"))

th3 <- theme_bw() +
  theme(text=element_text(size=30),
        legend.title=element_blank(),
        legend.text=element_text(size=10),
        legend.position="bottom",
        #text=element_text(family="Arial"),
        axis.title.x=element_text(margin=margin(10,0,0,0)),
        axis.title.y=element_text(margin=margin(0,10,0,0)),
        plot.margin=margin(1,1,1,1,"cm")) 
