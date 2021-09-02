library(ckbplotr) #it has to loaded at the beginning of the session
library("readxl")
library(data.table)
# http://ldsc.broadinstitute.org/lookup/
# Then I clicked on Genetic correlation results
# Then download the existing genetic correlation results for 221 traits (without 7 traits from ENIGMA) using data from LD Hub,
data <- read_excel("/home/gagelo01/workspace/Projects/Nooshin_NAFLD/Data/Raw/LD-Hub_genetic_correlation_221x221_no_ENIGMA.xlsx",
                    sheet = 1)

data<-data[,-1]
data <- as.matrix(data)

create_full_mat <- function(data,index){
new <-data
new[]<-NA
new[lower.tri(new)]<- map_chr(strsplit(data[lower.tri(data)], "\\s+"), `[`,index) #rg
new[upper.tri(new)] <- t(new)[upper.tri(new)] #this is the brilliant part
diag(new) <- 1
new <- as.data.frame(new)
new <- sapply(new, as.numeric)
new <- as.data.frame(new)
return(new)
}

undebug(create_full_mat)
data_rg <- create_full_mat(data, 3) #rg
data_se <- create_full_mat(data, 4) #se
data_p <- create_full_mat(data, 6) #p

index<- grep("ibd", tolower(colnames(data)))

data_IBD_full <-data.frame(traits = colnames(data_rg),est = data_rg[,index], se = data_se[,index], p = data_p[,index])


#create a forest plot of significant
setDT(data_IBD_full)
data_IBD <- data_IBD_full[p<0.05,]
data_IBD<- data_IBD[order(-est),]


resultsA <- data.frame(variable = substr(data_IBD$traits, 1, nchar(data_IBD$traits) -12),
                       estimate = data_IBD$est,
                       stderr =  data_IBD$se)
make_forest_plot(panels = list(resultsA),
                 col.key        = "variable",
                 exponentiate =  FALSE,
                 panel.headings = paste0("Genetic Correlation with IBD and ", nrow(data_IBD_full), " traits from LDHub (only p<0.05 shown)"),
                 nullval = 0,
                 col.right.heading = "Rg (95% CI)")

ggsave("/home/gagelo01/workspace/Projects/Nooshin_NAFLD/Results/Forest_plot_IBD_Rg.png",
       width=12,height=12,units="in",scale=1,
       device = "png")
