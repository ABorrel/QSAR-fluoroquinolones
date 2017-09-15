library(ggplot2)



multipleHist = function(d, prout){


  ggplot(d, aes(x=Pseudomonas.aeruginosa)) + 
      geom_histogram(aes(y=..density..),     # Histogram with density instead of count on y-axis
                     binwidth=.5,
                     colour="black", fill="white") +
    labs(x = "-log10(MIC)", y = "Frequencies") + 
  xlim (c(-2.5, 2.5))+
    theme(axis.text.y = element_text(size = 15, hjust = 0.5, vjust =0.1), axis.text.x = element_text(size = 15, hjust = 0.5, vjust =0.1), axis.title.y = element_text(size = 15, hjust = 0.5, vjust =0.1), axis.title.x =  element_text(size = 15, hjust = 0.5, vjust =0.1))+
      geom_density(alpha=.2, fill="#FF6666")
  ggsave(paste(prout, "Pseudomonas-aeruginosa.png", sep = ""), dpi = 300, width = 8, height = 7)
  
  ggplot(d, aes(x=Staphylococcus.aureus)) + 
      geom_histogram(aes(y=..density..),     # Histogram with density instead of count on y-axis
                     binwidth=.5, 
                     colour="black", fill="white") +
    labs(x = "-log10(MIC)", y = "Frequencies") + 
  xlim (c(-2.5, 2.5))+
    theme(axis.text.y = element_text(size = 15, hjust = 0.5, vjust =0.1), axis.text.x = element_text(size = 15, hjust = 0.5, vjust =0.1), axis.title.y = element_text(size = 15, hjust = 0.5, vjust =0.1), axis.title.x =  element_text(size = 15, hjust = 0.5, vjust =0.1))+
      geom_density(alpha=.2, fill="#A5FFB3")
  ggsave(paste(prout, "Staphyococcus-aureus.png", sep = ""), dpi = 300, width = 8, height = 7)
  
  ggplot(d, aes(x=Streptococcus.pneumoniae)) + 
      geom_histogram(aes(y=..density..),     # Histogram with density instead of count on y-axis
                     binwidth=.5,
                     colour="black", fill="white") +
    labs(x = "-log10(MIC)", y = "Frequencies") + 
  xlim (c(-2.5, 2.5))+
    theme(axis.text.y = element_text(size = 15, hjust = 0.5, vjust =0.1), axis.text.x = element_text(size = 15, hjust = 0.5, vjust =0.1), axis.title.y = element_text(size = 15, hjust = 0.5, vjust =0.1), axis.title.x =  element_text(size = 15, hjust = 0.5, vjust =0.1))+
      geom_density(alpha=.2, fill="#FAFFA5")
  ggsave(paste(prout, "Streptococcus-pneumoniae.png", sep = ""), dpi = 300, width = 8, height = 7)
  
  ggplot(d, aes(x=Escherichia.coli)) + 
      geom_histogram(aes(y=..density..),     # Histogram with density instead of count on y-axis
                     binwidth=.5,
                     colour="black", fill="white") +
    labs(x = "-log10(MIC)", y = "Frequencies") + 
  xlim (c(-2.5, 2.5))+
    theme(axis.text.y = element_text(size = 15, hjust = 0.5, vjust =0.1), axis.text.x = element_text(size = 15, hjust = 0.5, vjust =0.1), axis.title.y = element_text(size = 15, hjust = 0.5, vjust =0.1), axis.title.x =  element_text(size = 15, hjust = 0.5, vjust =0.1))+
      geom_density(alpha=.2, fill="#BAA5FF")
  ggsave(paste(prout, "Escherichia-coli.png", sep = ""), dpi = 300, width = 8, height = 7)
}