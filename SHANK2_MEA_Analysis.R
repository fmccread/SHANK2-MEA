setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) 
bd = getwd()

library(tibble)
library(plyr)
library(readr)
library(tidyr)
library(dplyr)
library(data.table)
library(ggplot2)
library(ggExtra)
library(stringr)
library(plotly)
library(kSamples)
library(pals)
library(cowplot)
library(cluster)
library(scatterplot3d)
library(parallel)
library(moments)
library(scico)
library(ggfortify)
source("MEA_metric_toolbox.R")



###########################
#         SETUP           #
###########################

#load plate metadata
shank2.info <- read.csv("SHANK2_recording_info.csv")

#comparisons for computing stats
line.comparisons = list(c("CTRL", "KO"), 
                        c("R841X-C", "R841X"), 
                        c("R841X", "R841X DHPG"), 
                        c("R841X-C", "R841X DHPG"), 
                        c("R841X-C", "R841X-C DHPG"))

# groups for plotting metrics
plot.groups1 = list(c("CTRL", "KO"), 
                    c("R841X-C", "R841X"),
                    c("R841X-C", "R841X", "R841X DHPG")) %>% setNames(., c("KO", "R841X", "DHPG"))

plot.groups2 = list(c("R841X", "R841X DHPG"),
                    c("R841X-C","R841X-C DHPG"),
                    c("R841X-C", "R841X", "R841X DHPG"),
                    c("R841X-C", "R841X-C DHPG", "R841X", "R841X DHPG")) %>% setNames(., c("DHPG - R841X", "DHPG - R841X-C",  "DHPG - Three", "DHPG - All"))

plot.groups3 = list(c("CTRL", "KO"), 
                    c("R841X-C", "R841X"),
                    c("R841X-C", "R841X-C DHPG", "R841X", "R841X DHPG")) %>% setNames(., c("KO", "R841X", "DHPG - All"))



# groups for plotting ISI distributions
isi.groups =  list(c("CTRL", "KO"), 
                   c("R841X-C", "R841X"),
                   c("R841X-C","R841X-C DHPG"),
                   c("R841X-C", "R841X DHPG")) %>% setNames(., c("KO", "R841X", "R841X DHPG", "R841X-C DHPG"))


plot.groups = list(c("R841X-C", "R841X DHPG")) %>% setNames(., c("Treated vs Untreated"))




################################
#         DATA CLEANUP         #
################################

#compile metrics
full.metrics <- read.csv("SHANK2_NeuralMetrics.csv")

#tidy and get averages
metric.averages = get_averages(full.metrics)



###########################
#         STATS           #
###########################

wilcox.stats = compile_stats_wilcox(tidy.metrics, line.comparisons, weeks = c(4,5,6,7,8))



#####################################
#         PLOTTING METRICS          #
#####################################

#plotting simple line plots
plot_simple_lines(full.metrics, groups = c("R841X", "R841X-C", "CTRL", "KO"))

#plotting pooled data
plot_pooled_lines(metric.averages,
                  groups = plot.groups1,
                  weeks = c(4,5,6,7,8),
                  stats = AD.stats,
                  #error.type = "CI",
                  error.type = "SEM",
                  #error.viz = "errorbar",
                  error.viz = "shade",
                  width = 6,
                  height = 5
)
plot_pooled_lines(metric.averages,
                  groups = plot.groups2,
                  weeks = c(4,5,6,7,8),
                  stats = wilcox.stats,
                  #error.type = "CI",
                  error.type = "SEM",
                  #error.viz = "errorbar",
                  error.viz = "shade",
                  width = 6.5,
                  height = 3.5
)

#plotting data by week
plot_weekly_dotplots(full.metrics, 
                     stats = wilcox.stats, 
                     groups = plot.groups3, 
                     weeks = c(7), 
                     width = 3.5, 
                     height = 5)


plot_weekly_ecdf(tidy.metrics, 
                 groups = plot.groups3, 
                 weeks = c(7),
                 width = 4,
                 height = 6)





#################################################
#         CorSE ANALYSIS AND PLOTTING           #
#################################################

#load data
corfiles = list.files(path = paste(bd, "/CorSE_Analysis", sep = ""), pattern = "CorrData.csv", recursive=TRUE)

#get file data
corse.info = data.frame()
for (i in seq_along(corfiles)) {
  
  #get file
  corfile = corfiles[i]
  
  #get metadata
  file.ID = str_split(corfile, "/") %>% unlist(.) %>% tail(.,1)
  plate = str_split(file.ID, "_") %>% unlist(.) %>% .[1] %>% gsub("Plate-", "", .) %>% as.numeric(.)
  day = str_split(file.ID, "_") %>% unlist(.) %>% .[2] %>% gsub("D", "", .) %>% as.numeric(.)
  week = day2week_v2(day)
  well = str_split(file.ID, "_") %>% unlist(.) %>% .[3]
  line = str_split(file.ID, "_") %>% unlist(.) %>% .[4]
  
  meta.df = data.frame(Plate = plate,
                       Day = day,
                       Week = week,
                       Well = well,
                       Line = line,
                       File.ID = file.ID,
                       Path = paste(bd, "/CorSE_Analysis/",corfile, sep = "")
  )
  corse.info = rbind(meta.df, corse.info)
  
}

#trim file data
corse.keep = filter(corse.info, Week %in% c(4,5,6,7,8))
corse.keep = filter(corse.keep, Line %in% c("R841X", "R841X-C", "CTRL", "KO", "R841X DHPG"))
corse.keep = filter(corse.keep, Plate %in% shank2.info$Plate)

#threshold for calling 'strong' positive correlations
thresh = 0.5

#getting densities of synchrony scores for each well
corse.densities = data.frame()
for (i in 1:nrow(corse.keep)) {
  
  cordata = read_csv(corse.keep$Path[i], col_names = FALSE, show_col_types = FALSE)
  melt.cordata = unlist(cordata) %>% unname(.)
  melt.cordata = melt.cordata[-c(which(melt.cordata == 0))]
  
  
  
  mean.cor = mean(melt.cordata)
  n.strong = length(melt.cordata[which(melt.cordata > thresh)])
  
  
  cor.density = density(melt.cordata, from = -0.5, to = 1, n = 700)
  
  p.df = data.frame(Synchrony.Score = cor.density$x,
                    Density = cor.density$y, 
                    Line = corse.keep$Line[i])
  
  ggplot(data = p.df, aes(x = Synchrony.Score, y = Density,  fill = Line)) +
    geom_area(color= "black")+
    geom_line(lwd = 1.25, color = "black") +
    scm + 
    sfm + 
    ylab("Average Density") +
    xlab("CorSE Synchrony Score")+
    theme(
      axis.title.y = element_text( size = 12,face = "bold", color = 'black'),
      axis.title.x = element_text( size = 12,face = "bold", color = 'black'),
      axis.text.y = element_text(size = 12, face = "bold",color = 'black'),
      axis.text.x = element_text(size = 12, face = "bold",color = 'black'),
      panel.background = element_rect(fill = 'transparent'),
      axis.line = element_line(size = 1, color = "black"),
      legend.title = element_blank(),
      legend.position = "none"
    )
  sv.path = gsub(corse.keep$File.ID[i], "", corse.keep$Path[i])
  sv.name = gsub("CorrData.csv", "Synchrony-Score-Density.pdf", corse.keep$File.ID[i])
  ggsave(paste(sv.path, "/", sv.name, sep = ""),height = 5, width = 7)
  
  temp.df = cbind(corse.keep[i,c(1:5)], data.frame(Mean.CorSE = mean.cor, n.Strong = n.strong), as.data.frame(t(cor.density$y)))
  corse.densities = rbind(corse.densities, temp.df)
  
}
z = cor.density$x

##plotting

lines = unique(corse.densities$Line)

#get average densities
avg.corse.densities = data.frame()

for (a in seq_along(lines)) {
  
  #subset data
  p.line = lines[a]
  corse.dat = filter(corse.densities, Line == p.line & Week %in% c(7,8))[,-c(1:7)]
  corse.avg = colMeans(corse.dat)
  
  
  avg.df = data.frame(Synchrony.Score = z, Density = corse.avg, Line = p.line)
  avg.corse.densities = rbind(avg.corse.densities, avg.df)
  
}

comps = list(c("CTRL", "KO"), c("R841X", "R841X-C"), c("R841X-C", "R841X DHPG"), c("R841X-C", "R841X-C DHPG"), c("R841X", "R841X DHPG"))

for (p in seq_along(comps)) {
  
  #subset data
  p.comp = comps[[p]]
  plot.df = filter(avg.corse.densities, Line %in% p.comp)
  plot.df$Line = factor(plot.df$Line, levels = c(p.comp[1], p.comp[2]))
  
  ggplot(data = plot.df, aes(x = Synchrony.Score, y = Density, color = Line, fill = Line)) +
    geom_line(lwd = 1.25) +
    scm + 
    sfm + 
    ylab("Average Density") +
    xlab("CorSE Synchrony Score")+
    ylim(0,5)+
    theme(
      axis.title.y = element_text( size = 12,face = "bold", color = 'black'),
      axis.title.x = element_text( size = 12,face = "bold", color = 'black'),
      axis.text.y = element_text(size = 12, face = "bold",color = 'black'),
      axis.text.x = element_text(size = 12, face = "bold",color = 'black'),
      panel.background = element_rect(fill = 'transparent'),
      axis.line = element_line(size = 1, color = "black"),
      legend.title = element_blank(),
      legend.position = "none"
    )
  ggsave(paste("Synchrony_Score_Density_", p, ".pdf", sep =""), height = 4, width = 4)
  
  
  ggplot(data = plot.df, aes(x = Synchrony.Score, y = Density, color = Line, fill = Line)) +
    geom_line(lwd = 1.25) +
    scm + 
    sfm + 
    ylab("Average Density") +
    xlab("CorSE Synchrony Score")+
    ylim(0,0.6)+
    xlim(0.5,1)+
    theme(
      axis.title.y = element_text( size = 12,face = "bold", color = 'black'),
      axis.title.x = element_text( size = 12,face = "bold", color = 'black'),
      axis.text.y = element_text(size = 12, face = "bold",color = 'black'),
      axis.text.x = element_text(size = 12, face = "bold",color = 'black'),
      panel.background = element_rect(fill = 'transparent'),
      axis.line = element_line(size = 1, color = "black"),
      legend.title = element_blank(),
      legend.position = "none"
    )
  ggsave(paste("Synchrony_Score_Density_expanded_", p, ".pdf", sep =""), height = 5, width = 4)
  
}

