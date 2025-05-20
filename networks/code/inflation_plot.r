library(ggplot2, lib.loc="/Storage/data1/jorge.munoz/NRGSC.new/libraries")

df <- read.table("/Storage/data1/jorge.munoz/NRGSC.new/networks/results/I_test.csv", header = FALSE, sep = ",")
colnames(df) <- c("Inflation", "No_clusters")

plot <- ggplot(df, aes( x = Inflation, y = No_clusters)) +
  geom_point(lwd=2.6, colour = "black")+
  geom_line(lwd=1)+
  labs(title="Inflation value evaluation",x="Inflation", y = "Number of clusters")+
  theme_bw() +
  theme( text = element_text(family = "Times New Roman", size=20),
         legend.key.size = unit(1, 'cm'), 
         panel.border = element_blank(),
         panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         axis.text.x=element_text(colour="black" , hjust=0.5),
         axis.text.y=element_text(colour="black"),
         axis.line = element_line(colour = "black", size = 1.8))

ggsave("/Storage/data1/jorge.munoz/NRGSC.new/networks/results/inflation_test.png", plot, device = "png",  units = "cm", width =  18, height = 12)

