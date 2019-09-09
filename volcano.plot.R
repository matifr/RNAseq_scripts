volcano.plot = function(tt){
  library(ggplot2)
  tt$threshold = as.factor(abs(tt$logFC) >= 1.5 & tt$P.Value < 0.05)
  g = ggplot(data= as.data.frame(tt), aes(x=logFC, y=-log10(P.Value), colour=threshold)) +
    geom_point(alpha=0.5, size=1.5)+
    theme(axis.title.x=element_text(size=25), legend.position = "none",
          axis.title.y=element_text(size=25)) +
    theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
    theme(plot.title = element_text(face="bold", size=20)) +
    theme(axis.text=element_text(size=12)) +
    theme(text = element_text(size=10)) +
    theme(axis.line = element_line(size = 0.5, colour = "black"))+
    theme(panel.border = element_rect(colour = "white",size = 0.1))+
    xlim(c(-7,7)) +
    
    xlab("log2 fold change")  +
    ylab("-log10 p-value")
  g
}