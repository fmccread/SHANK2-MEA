
##### General Helper Functions #####
day2week <- function(x){
  
  #Convert days to weeks
  day.in = seq(1,81,by=1)
  day.out = c(rep(0,4), rep(0.5,3), rep(1,4), rep(1.5,3), rep(2,4), rep(2.5,3), rep(3,4), rep(3.5,3), rep(4,4), rep(4.5,3), rep(5,4), rep(5.5,3),
              rep(6,4), rep(6.5,3), rep(7,4), rep(7.5,3), rep(8,4), rep(8.5,3), rep(9,4), rep(9.5,3), rep(10,4), rep(10.5,3), rep(11,4))
  
  x <- plyr::mapvalues(x, from = day.in, to = day.out, warn_missing = FALSE)
}

day2week_v2 <- function(x){
  
  #Convert days to weeks
  day.in = seq(1,84,by=1)
  day.out = c(rep(0,7), rep(1,7), rep(2,7), rep(3,7), rep(4,7), rep(5,7), rep(6,7), rep(7,7), rep(8,7), rep(9,7), rep(10,7), rep(11,7))
  
  x <- plyr::mapvalues(x, from = day.in, to = day.out, warn_missing = FALSE)
}

se <- function(x) {
  sqrt(var(x[!is.na(x)])/length(x[!is.na(x)]))
}

`%notin%` <- Negate(`%in%`)




##### For Compiling Neural Metrics #####
split_neuralMetrics <-function(bd){
  
  
  files = list.files(path = paste(bd,"Spike Data", sep = "/"), pattern = "neuralMetrics", recursive = TRUE)
  dir.create(paste(bd,"Neural Metrics",sep="/"), showWarnings = FALSE)
  
  
  for (i in seq_along(files)) {
    
    file = paste("Spike Data/", files[i], sep = "")
    metrics = read_neuralMetrics(file)
    
    file = gsub("\\(6 x STD\\)_\\(000\\)", "\\(6 x STD\\)\\(000\\)", file)
    metadat = get_metadata(file)
    
    plate = metadat[[1]]
    day = metadat[[2]]
    week = metadat[[3]]
    algo = metadat[[4]]
    
    dir.create(paste("Neural Metrics/Plate ", plate, sep = ""), showWarnings = FALSE)
    
    wells = colnames(metrics)[-1]
    
    
    for (j in seq_along(wells)) {
      
      #subset well
      well = wells[j]
      sub.df = select(metrics, Metric, well)
      
      #save
      fname = paste("Plate-", plate, "_D", day, "_", well, "_NeuralMetrics_", algo, ".csv", sep = "")
      fpath = paste("Neural Metrics/Plate ", plate, "/", fname, sep = "")
      write.csv(sub.df, file = fpath, row.names = FALSE)
      
    }
    
  }
  
}

compile_neuralMetrics <- function(bd, shank2.info){
  
  #split neural metrics files
  #message("splitting NeuralMetrics files...")
  #split_neuralMetrics(bd)
  
  #compile neural metrics files
  message("Compiling optimized NeuralMetrics files...")
  well.keep = filter(shank2.info, Analyze == "YES")
  full.metrics = data.frame()
  
  for (i in 1:nrow(well.keep)) {
    
    #extract metadata
    plate = well.keep$Plate[i]
    day = well.keep$Day[i]
    well = well.keep$Well[i]
    line = well.keep$Line[i]
    algo = well.keep$Algorithm[i]
    week = day2week(day)
    
    #create metadata df
    meta.df = data.frame(Plate = plate,
                         Day = day,
                         Week = week,
                         Well = well,
                         Line = line,
                         NB.Algorithm = algo)
    
    #neuralmetric file path from metadata
    fname = paste("Plate-", plate, "_D", day, "_", well, "_NeuralMetrics_", algo, ".csv", sep = "")
    fpath = paste(bd, "/Neural Metrics/", "Plate ", plate, "/", fname, sep = "")
    
    #load file
    data = read_csv(fpath, show_col_types = FALSE, progress = FALSE)[-1,] %>% t(.) %>% as.data.frame(.)
    
    #clean up data
    colnames(data) = data[1,] %>% make.names(.) %>% clean_names(.)
    data = data[-1,]
    data = cbind(meta.df, data)
    
    #append
    full.metrics = rbind(full.metrics, data)
    
  }
  
  #fix line names  
  full.metrics$Line = gsub("19-2-2 WT", "CTRL", full.metrics$Line)
  full.metrics$Line = gsub("19-2-2 KO", "KO", full.metrics$Line)
  
  #save data
  write.csv(full.metrics, paste(bd, "/Neural Metrics/SHANK2_Compiled_NeuralMetrics.csv", sep = ""), row.names = FALSE)
  
  return(full.metrics)
  
}

compile_neuralMetrics_ISI <- function(bd, shank2.info){
  
  #split neural metrics files
  #message("splitting NeuralMetrics files...")
  #split_neuralMetrics(bd)
  
  #compile neural metrics files
  message("Compiling optimized NeuralMetrics files... ISI")
  well.keep = filter(shank2.info, Analyze == "YES")
  full.metrics = data.frame()
  
  for (i in 1:nrow(well.keep)) {
    
    #extract metadata
    plate = well.keep$Plate[i]
    day = well.keep$Day[i]
    well = well.keep$Well[i]
    line = well.keep$Line[i]
    algo = "ISI"
    week = day2week(day)
    ID = paste(line, algo, sep = "_")
    
    #create metadata df
    meta.df = data.frame(Plate = plate,
                         Day = day,
                         Week = week,
                         Well = well,
                         Line = line,
                         NB.Algorithm = algo,
                         ID = ID)
    
    #neuralmetric file path from metadata
    fname = paste("Plate-", plate, "_D", day, "_", well, "_NeuralMetrics_ISI.csv", sep = "")
    fpath = paste(bd, "/Neural Metrics/", "Plate ", plate, "/", fname, sep = "")
    
    #load file
    data = read_csv(fpath, show_col_types = FALSE, progress = FALSE)[-1,] %>% t(.) %>% as.data.frame(.)
    
    #clean up data
    colnames(data) = data[1,] %>% make.names(.) %>% clean_names(.)
    data = data[-1,]
    data = cbind(meta.df, data)
    
    #append
    full.metrics = rbind(full.metrics, data)
    
  }

  #fix line names  
  full.metrics$Line = gsub("19-2-2 WT", "CTRL", full.metrics$Line)
  full.metrics$Line = gsub("19-2-2 KO", "KO", full.metrics$Line)
  
  #save data
  write.csv(full.metrics, paste(bd, "/Neural Metrics/SHANK2_Compiled_NeuralMetrics_ISI.csv", sep = ""), row.names = FALSE)
  
  return(full.metrics)
  
}

compile_neuralMetrics_Envelope <- function(bd, shank2.info){
  
  #split neural metrics files
  #message("splitting NeuralMetrics files...")
  #split_neuralMetrics(bd)
  
  #compile neural metrics files
  message("Compiling optimized NeuralMetrics files... Envelope")
  well.keep = filter(shank2.info, Analyze == "YES")
  full.metrics = data.frame()
  
  for (i in 1:nrow(well.keep)) {
    
    #extract metadata
    plate = well.keep$Plate[i]
    day = well.keep$Day[i]
    well = well.keep$Well[i]
    line = well.keep$Line[i]
    algo = "Envelope"
    week = day2week(day)
    ID = paste(line, algo, sep = "_")
    
    #create metadata df
    meta.df = data.frame(Plate = plate,
                         Day = day,
                         Week = week,
                         Well = well,
                         Line = line,
                         NB.Algorithm = algo,
                         ID = ID)
    
    #neuralmetric file path from metadata
    fname = paste("Plate-", plate, "_D", day, "_", well, "_NeuralMetrics_Envelope.csv", sep = "")
    fpath = paste(bd, "/Neural Metrics/", "Plate ", plate, "/", fname, sep = "")
    
    #load file
    data = read_csv(fpath, show_col_types = FALSE, progress = FALSE)[-1,] %>% t(.) %>% as.data.frame(.)
    
    #clean up data
    colnames(data) = data[1,] %>% make.names(.) %>% clean_names(.)
    data = data[-1,]
    data = cbind(meta.df, data)
    
    #append
    full.metrics = rbind(full.metrics, data)
    
  }
  
  #fix line names  
  full.metrics$Line = gsub("19-2-2 WT", "CTRL", full.metrics$Line)
  full.metrics$Line = gsub("19-2-2 KO", "KO", full.metrics$Line)
  
  #save data
  write.csv(full.metrics, paste(bd, "/Neural Metrics/SHANK2_Compiled_NeuralMetrics_Envelope.csv", sep = ""), row.names = FALSE)
  
  return(full.metrics)
  
}

compile_neuralMetrics_Adaptive <- function(bd, shank2.info){
  
  #split neural metrics files
  #message("splitting NeuralMetrics files...")
  #split_neuralMetrics(bd)
  
  #compile neural metrics files
  message("Compiling optimized NeuralMetrics files... Adaptive")
  well.keep = filter(shank2.info, Analyze == "YES")
  full.metrics = data.frame()
  
  for (i in 1:nrow(well.keep)) {
    
    #extract metadata
    plate = well.keep$Plate[i]
    day = well.keep$Day[i]
    well = well.keep$Well[i]
    line = well.keep$Line[i]
    algo = "Adaptive"
    week = day2week(day)
    ID = paste(line, algo, sep = "_")
    
    #create metadata df
    meta.df = data.frame(Plate = plate,
                         Day = day,
                         Week = week,
                         Well = well,
                         Line = line,
                         NB.Algorithm = algo,
                         ID = ID)
    
    #neuralmetric file path from metadata
    fname = paste("Plate-", plate, "_D", day, "_", well, "_NeuralMetrics_adaptive.csv", sep = "")
    fpath = paste(bd, "/Neural Metrics/", "Plate ", plate, "/", fname, sep = "")
    
    #load file
    data = read_csv(fpath, show_col_types = FALSE, progress = FALSE)[-1,] %>% t(.) %>% as.data.frame(.)
    
    #clean up data
    colnames(data) = data[1,] %>% make.names(.) %>% clean_names(.)
    data = data[-1,]
    data = cbind(meta.df, data)
    
    #append
    full.metrics = rbind(full.metrics, data)
    
  }
  
  #fix line names  
  full.metrics$Line = gsub("19-2-2 WT", "CTRL", full.metrics$Line)
  full.metrics$Line = gsub("19-2-2 KO", "KO", full.metrics$Line)
  
  #save data
  write.csv(full.metrics, paste(bd, "/Neural Metrics/SHANK2_Compiled_NeuralMetrics_Adaptive.csv", sep = ""), row.names = FALSE)
  
  return(full.metrics)
  
}

read_neuralMetrics <- function(file){
  
  #for loading file
  if(str_detect(file, "Adaptive")){
    skip = 25
  } else if (str_detect(file, "Envelope")){
    skip = 26
  } else if (str_detect(file, "ISI-200")){
    skip = 25
  } else if (str_detect(file, "adaptive")){
    skip = 25
  } else if (str_detect(file, "ISI.csv")){
    skip = 26
  }
  
  #assigning any names to be changed/standardized
  names.in = c("ex.7 h6", "ex.7 H6", "Ex.7-H6", "EX.7 H6 (-/-)",
               "SK19-2-2", "SK0019-2-2", "SK0019-002 #2",
               "Corr B8-6A", "CORR B8-6A",
               "SK441-3-5", "SK0441-003 #5",
               "CORR B8-6A DHPG [10 µM]", "CORR B8-6A DHPG [10 μM]",
               "SK0441-003 #5 DHPG [10 µM]", "SK0441-003 #5 DHPG [10 μM]",
               "CORR B8-6A DHPG [20 µM]",
               "SK0441-003 #5 DHPG [20 µM]",
               "Treatment/ID"
  )
  
  names.out = c(rep("19-2-2 KO", 4), 
                rep("19-2-2 WT",3),
                rep("R841X-C", 2),
                rep("R841X", 2),
                rep("R841X-C_DHPG",2),
                rep("R841X_DHPG", 2),
                rep("R841X-C_DHPG", 1),
                rep("R841X_DHPG", 1),
                "Cell_Line")
  
  
  #load file
  metrics <- read.csv(file, skip = skip)
  
  #cleanup
  metrics$Well.Averages = gsub("    ", "", metrics$Well.Averages)
  metrics = metrics[1:(which(metrics$Well.Averages == "Measurement")-1),]
  metrics = metrics[-c(which(metrics$Well.Averages == "Activity Metrics")),]
  metrics = metrics[-c(which(metrics$Well.Averages == "Electrode Burst Metrics")),]
  metrics = metrics[-c(which(metrics$Well.Averages == "Network Burst Metrics")),]
  metrics = metrics[-c(which(metrics$Well.Averages == "Synchrony Metrics")),]
  metrics = metrics[-c(which(metrics$Well.Averages == "Average Network Burst Metrics")),]
  metrics = metrics[-c(which(metrics$Well.Averages == "Network Bursts Ignored Flag")),]
  metrics = metrics[,-c(which(colnames(metrics) == "X"))]
  
  #fix line names
  lines = metrics[1,] %>% unlist()
  lines = plyr::mapvalues(lines, from = names.in, to = names.out, warn_missing = FALSE)
  metrics[1,] = lines
  
  #fix column names
  colnames(metrics)[1] <- "Metric"
  
  return(metrics)
  
}

get_metadata <- function(file){
  
  split.fname = str_split(file, "_") %>% unlist()
  
  plate = split.fname[max(which(str_detect(split.fname, "Plate")))] %>% gsub("Plate", "",.) %>% gsub("-", "",.) %>% as.numeric(.)
  day = split.fname[which(str_detect(split.fname, "D\\d\\d"))] %>% gsub("D", "",.) %>% as.numeric(.)
  week = day2week(day)
  algo = split.fname[which(str_detect(split.fname, ".csv"))] %>% gsub(".csv", "",.)
  
  metadat = list(plate, day, week, algo)
  
  return(metadat)
  
}

clean_names <- function(m.names){
  
  m.names = gsub("\\.\\.\\.", "\\.", m.names)
  m.names = gsub("\\.\\.", "\\.", m.names)
  m.names = gsub("Hz\\.", "Hz", m.names)
  m.names = gsub("sec\\.", "sec", m.names)
  m.names = gsub("ms\\.", "ms", m.names)
  
  return(m.names)
  
}

tidy_neuralMetrics <- function(full.metrics){
  
  #collapse weeks
  full.metrics$Week <- floor(full.metrics$Week)

  #correct names
  full.metrics$Line = gsub("19-2-2 WT", "CTRL", full.metrics$Line)
  full.metrics$Line = gsub("19-2-2 KO", "KO", full.metrics$Line)
  
  #plates to keep
  #ko.plates = c(14,16,36,38)
  ko.plates = c(11, 14, 16, 24, 36, 38, 41)
  mut.plates= c(31,32,35,37,39,47,48,49,50)
  dhpg.plates = c(32,35,43,45,47,48,49,50)
  plates = c(ko.plates, mut.plates, dhpg.plates) %>% unique(.)
  
  #filter
  trimmed.metrics = full.metrics[,-c(2,4,6,7,11,13,14,17,18,21,22,26,27,28,31,32,36,38,40,41,45,46,49,50,54,55,56,58,59,60,61,63,64,65,70,72,75,76)]
  trimmed.metrics = filter(trimmed.metrics, Plate %in% plates)
  #trimmed.metrics = filter(trimmed.metrics, Plate %in% c(31,32,35,37,39,47,48,49,50,40,43,45,14,16,36,38,12,15,11,24,30,33,41))
  
  #trimmed.metrics = trimmed.metrics[-c(which(trimmed.metrics$Line %in% c("R841X DHPG", "R841X-C DHPG") & !(trimmed.metrics$Plate %in% c(40,43,45)))),]
  
  return(trimmed.metrics)
  
}

export_classifier_metrics <- function(tidy.metrics){
  
  #subset weeks to plot
  tidy.metrics = filter(tidy.metrics, Week %in% c(4,5,6,7,8))
  
  #subset lines
  R841X.ml = filter(tidy.metrics, Line %in% c("R841X", "R841X-C"))
  DHPG.ml = filter(tidy.metrics, Line == "R841X DHPG")
  ko.ml = filter(tidy.metrics, Line %in% c("CTRL", "KO"))
  
  #save
  write.csv(R841X.ml, "Classifier Metrics/neuralMetrics_Mut-vs-Cor.csv", row.names = FALSE)
  write.csv(ko.ml, "Classifier Metrics/neuralMetrics_KO-vs-WT.csv", row.names = FALSE)
  write.csv(DHPG.ml, "Classifier Metrics/neuralMetrics_DHPG.csv", row.names = FALSE)
  
}



##### For Computing Stats #####

get_averages <- function(tidy.metrics){
  
  metric.averages = data.frame()
  lines = unique(tidy.metrics$Line)
  
  for (i in seq_along(lines)) {
    
    #subset line
    line = lines[i]
    line.df = filter(tidy.metrics, Line == line)
    weeks = unique(line.df$Week)
    
    for (j in seq_along(weeks)) {
      
      #subset week
      week = weeks[j]
      week.df = filter(line.df, Week == week)
      mets = colnames(week.df)[-c(1:3)]
      
      for (k in seq_along(mets)) {
        
        #subset metric
        metric = mets[k]
        met.df = select(week.df, all_of(metric))
        met.df = na.omit(met.df)
        colnames(met.df) <- "Metric"
        
        #skip if all NA
        if(nrow(met.df) == 0){
          next
        }
        
        #get mean and SE
        met.mean = mean(met.df$Metric)
        met.se = se(met.df$Metric)
        se.yneg = met.mean - met.se
        se.ypos = met.mean + met.se
        
        #get 95% confidence interval
        lmod <- lm(Metric ~ 1, met.df)
        ci <- confint(lmod, level=0.95)
        ci.yneg = ci[1]
        ci.ypos = ci[2]
        
        #create data frame for results
        temp.df = data.frame(Line = line,
                             Week = week,
                             Metric = metric,
                             Value.Mean = met.mean,
                             Value.SE = met.se,
                             SE.yneg = se.yneg,
                             SE.ypos = se.ypos,
                             CI.yneg = ci.yneg,
                             CI.ypos = ci.ypos)
        
        #append
        metric.averages = rbind(metric.averages, temp.df)
        
      }
    }
  }
  
  write.csv(metric.averages, "NeuralMetrics_Averages_CI.csv")
  return(metric.averages)
}

compile_stats_wilcox <- function(full.metics, comparisons, weeks){
  
  mets = colnames(full.metics)[-c(1:3)]
  metric.df = select(full.metics, Week, Line, mets)
  wilcox.df = data.frame()
  
  #n for testing correction
  n = length(mets)
  
  for (i in seq_along(comparisons)) {
    
    lines = comparisons[[i]]
    ctrl.line = lines[1]
    comp.line = lines[2]
    
    group.df = filter(full.metics, Line %in% lines)
    
    for (j in seq_along(weeks)) {
      
      #subset week
      week = weeks[j]
      week.df = filter(group.df, Week == week)
      
      WC.week = data.frame()
      
      for (k in seq_along(mets)){
        
        #subset metric
        met = mets[k]
        met.lab = gsub("\\.", " ", met)
        metric.df =  select(week.df, Week, Line, met)
      
        #subset control group data
        ctrl.dat = filter(metric.df, Line == ctrl.line)[[3]] %>% .[!is.na(.)]
        ctrl.n = length(ctrl.dat)
        ctrl.mean = mean(ctrl.dat)
        ctrl.se = sd(ctrl.dat) / sqrt(ctrl.n)
        
        #subset comparison group data
        comp.dat = filter(metric.df, Line == comp.line)[[3]] %>% .[!is.na(.)]
        comp.n = length(comp.dat)
        comp.mean = mean(comp.dat)
        comp.se = sd(comp.dat) / sqrt(comp.n)
        
        #Enter NS if no data
        if (ctrl.n == 0 | comp.n == 0){
          
          temp.df = data.frame(Week = week,
                               Metric = met.lab,
                               Control.Group = ctrl.line,
                               Comparison.Group = comp.line,
                               Control.n = ctrl.n,
                               Comparison.n = comp.n,
                               Control.Mean = NA,
                               Control.SEM = NA,
                               Comparison.Mean = NA,
                               Comparison.SEM = NA,
                               Fold.Change = NA,
                               y.Offset = NA,
                               W.Statistic = NA,
                               p.value = NA,
                               Significance = NA,
                               p.adj.bon = NA,
                               Significance.Bon = NA)
          
        } else {
          
          #fold change and Wilcox stats
          fc = comp.mean / ctrl.mean
          wc = wilcox.test(ctrl.dat, comp.dat)
          w.stat = wc$statistic
          pval = wc$p.value
          
          #get significance star
          if (pval < 0.005){
            sig = "***"
          } else if (pval < 0.01){
            sig = "**"
          } else if (pval < 0.05){
            sig = "*"
          } else {
            sig = ""
          }
          
          
          #bonferroni correction
          p.adj.bon = p.adjust(pval, method = "bonferroni", n)
          
          #get significance star
          if (p.adj.bon < 0.005){
            sig.bon = "***"
          } else if (p.adj.bon < 0.01){
            sig.bon = "**"
          } else if (p.adj.bon < 0.05){
            sig.bon = "*"
          } else {
            sig.bon = ""
          }
          
          
          y.offset = max(c(ctrl.dat, comp.dat)) + max(c(ctrl.dat, comp.dat)) * 0.2
          
          #Create data frame
          temp.df = data.frame(Week = week,
                               Metric = met.lab,
                               Control.Group = ctrl.line,
                               Comparison.Group = comp.line,
                               Control.n = ctrl.n,
                               Comparison.n = comp.n,
                               Control.Mean = ctrl.mean,
                               Control.SEM = ctrl.se,
                               Comparison.Mean = comp.mean,
                               Comparison.SEM = comp.se,
                               Fold.Change = fc,
                               y.Offset = y.offset,
                               W.Statistic = w.stat,
                               p.value = pval,
                               Significance = sig,
                               p.adj.bon = p.adj.bon,
                               Significance.Bon = sig.bon)
          
        }
        #append
        WC.week = rbind(WC.week, temp.df)
      }
      #Apply BH correction by week
      WC.week$p.adj.BH = p.adjust(WC.week$p.value, method = "BH")
      WC.week$Significance.BH <- ""
      WC.week$Significance.BH[which(WC.week$p.adj.BH < 0.05)] <- "*"
      WC.week$Significance.BH[which(WC.week$p.adj.BH < 0.01)] <- "**"
      WC.week$Significance.BH[which(WC.week$p.adj.BH < 0.005)] <- "***"
      
      #append
      wilcox.df = rbind(wilcox.df, WC.week)
    }
  }
  write.csv(wilcox.df, "NeuralMetrics_Wilcox_stats.csv")
  return(wilcox.df)
}



##### colour scales #####

sfm <- scale_fill_manual("Genotype", breaks=c("R841X","R841X DHPG","R841X-C","R841X-C DHPG", "KO", "CTRL"), 
values = c("#7FC97F", "#36855d","#FDAE61", "#de754b","#386CB0", "#F0027F"))

scm <- scale_color_manual("Genotype", breaks=c("R841X","R841X DHPG","R841X-C","R841X-C DHPG", "KO", "CTRL"), 
values = c("#7FC97F", "#36855d","#FDAE61", "#de754b","#386CB0", "#F0027F"))




##### Plotting Functions #####

plot_simple_lines <- function(full.metrics, groups){
  
  #display message
  message("Plotting Simple Line Plots...")
  
  #collapse weeks
  full.metrics$Week <- floor(full.metrics$Week)

  #correct names
  full.metrics$Line = gsub("19-2-2 WT", "CTRL", full.metrics$Line)
  full.metrics$Line = gsub("19-2-2 KO", "KO", full.metrics$Line)
  
  #Clean up data
  full.metrics = filter(full.metrics, Plate %in% c(15,17,18,31,32,35,37,39,47,48,49,50,40,43,45,14,16,36,38))
  weeks = c(2,3,4,5,6,7,8)
  
  plot.metrics = colnames(full.metrics)[-c(1:6)]
  plot.data = select(full.metrics, Plate, Week, Line, plot.metrics) %>% filter(., Line %in% groups & Week %in% weeks)
  
  #calculate means for each group by week
  plot.means = data.frame()
  for (i in seq_along(weeks)) {
    
    week.df = filter(plot.data, Week == weeks[i])
    
    for (j in seq_along(groups)){
      
      group.df = filter(week.df, Line == groups[j])
      
      #calculate means
      means = transpose(as.data.frame(colMeans(na.omit(group.df[,-c(1:3)]))))
      colnames(means) <- colnames(group.df[,-c(1:3)])
      
      #calculate standard error
      sem = transpose(as.data.frame(sapply(na.omit(group.df[,-c(1:3)]),function(x)sd(x)/sqrt(length(x)))))
      colnames(sem) <- paste(colnames(group.df[,-c(1:3)]), "SE", sep = ".")
      
      #metadata
      meta.df = data.frame(Week = weeks[i],
                           Line = groups[j])
      
      #combine
      temp.df = cbind(meta.df, means, sem)
      plot.means = rbind(plot.means, temp.df)
      
    }
  }
  
  for (p in seq_along(plot.metrics)) {
    
    met = plot.metrics[p]
    met.se = paste(met, "SE", sep = ".")
    met.lab = gsub("\\.", " ", met)
    met.lab = gsub("Hz", "(Hz)", met.lab)
    met.lab = gsub("sec", "(sec)", met.lab)
    
    plot.df = select(plot.means, Week, Line, met, met.se)
    colnames(plot.df)[c(3:4)] <- c("Mean", "SE")
    
    ggplot(plot.df, aes(x = Week, y = Mean, color = Line)) +
      geom_line(lwd = 1.4)+
      geom_point(size = 2)+
      geom_errorbar(aes(ymin=Mean-SE, ymax=Mean+SE), lwd=1, width=0.15)+
      scm +
      ylab(met.lab)+
      xlab("Week")+
      theme(axis.title.y = element_text(size = 20, color = "black", face = "bold"),
            axis.title.x = element_text(size = 20, color = "black", face = "bold"),
            plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
            axis.text.y = element_text(size = 18, color = "black", face = "bold"),
            axis.text.x = element_text(size = 18, color = "black", face = "bold"),
            panel.background = element_blank(),
            panel.spacing = unit(0, "lines"),
            strip.background = element_blank(),
            strip.placement = "outside",
            strip.text = element_text(size = 20, color = "black", hjust = 0.5, face = "bold"),
            axis.line = element_line(color = "black", size =1),
            legend.title = element_blank(),
            legend.key = element_blank(),
            legend.text = element_text(color = "black", face = "bold", size = 12),
            panel.spacing.x = unit(2, "lines")
      )
    path.name = "Plotting Output/Line Plots - All Lines/"
    dir.create(path.name, showWarnings = FALSE)
    ggsave(paste(path.name, met.lab, ".pdf", sep = ""), width = 10, height = 5)
    
  }
  
}

plot_pooled_violins <- function(full.metics, wilcox.stats, groups, weeks){
  
  #display message
  message("Plotting Pooled Violin Plots...")
  
  #get names of metrics for plotting
  mets = colnames(full.metics)[-c(1:3)]
  
  #subset weeks to plot
  full.metics = filter(full.metics, Week %in% weeks)
  
  for (i in seq_along(groups)) {
    
    #subset group for comparison and plotting
    group = groups[[i]]
    
    for (j in seq_along(mets)) {
      
      #metric to plot
      met = mets[j]
      met.lab = gsub("\\.", " ", met)
      
      #filter data
      plot.df = filter(full.metics, Line %in% group) %>% select(., Plate, Week, Line, met)
      
      #filter stats
      if(length(group) == 2){
        plot.stats = filter(wilcox.stats, Control.Group == group[1] & Comparison.Group == group[2] & Metric == met.lab & Week %in% weeks)
      } else if(length(group) == 3){
        plot.stats = filter(wilcox.stats, Control.Group == group[1] & Comparison.Group %in% c(group[2], group[3]) & Metric == met.lab & Week %in% weeks)
      }
      
      #clean up data for plotting
      plot.df$Line = factor(plot.df$Line, levels = group)
      plot.df$Week = paste("Week", plot.df$Week)
      
      #clean up stats for plotting
      plot.stats$Line = plot.stats$Comparison.Group
      plot.stats$Week = paste("Week", plot.stats$Week)
      
      #clean up plot label
      met.lab = gsub("Avg", "", met.lab)
      met.lab = gsub("Number of", "", met.lab)
      met.lab = gsub("Within", "in", met.lab)
      met.lab = gsub("Cross Correlation", "CC", met.lab)
      met.lab = gsub("Full Width at Half Height of", "Width of", met.lab)
      
      #plot
      ggplot(plot.df, aes_string(y=met, x = "Line", fill = "Line")) + 
        geom_violin(trim = FALSE, size = 1, scale = "width", color = 'black') + 
        geom_text(data = plot.stats, aes(x = Line, y = y.Offset, label = Significance), size = 10) +
        facet_wrap(.~ Week, strip.position = "bottom", scales = "free_x", ncol = 5) +
        stat_summary(fun = mean, fun.min = mean, fun.max = mean, geom = "crossbar", color = "black", width = 1, lwd = 0.7) +
        sfm +
        scm +
        ylab(met.lab)+
        theme(axis.title.y = element_text(size = 20, color = "black", face = "bold"),
              axis.title.x = element_blank(),
              plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
              axis.text.y = element_text(size = 18, color = "black", face = "bold"),
              axis.text.x = element_blank(),
              axis.ticks.x = element_blank(),
              panel.background = element_blank(),
              panel.spacing = unit(0, "lines"),
              strip.background = element_blank(),
              strip.placement = "outside",
              strip.text = element_text(size = 20, color = "black", hjust = 0.5, face = "bold"),
              axis.line = element_line(color = "black", size =1),
              legend.position = "none",
              panel.spacing.x = unit(2, "lines")
        )
      
      #create directory and save
      path.name = "Plotting Output/Violin Plots - Pooled/"
      dir.create(path.name, showWarnings = FALSE)
      path.name = paste(path.name, names(groups[i]), sep = "")
      dir.create(path.name, showWarnings = FALSE)
      ggsave(paste(path.name, "/", met.lab, ".pdf", sep = ""), width = 10, height = 5)
      
    }
    
  }
  
}

plot_weekly_dotplots <- function(full.metics, stats, groups, weeks, width, height){
  
  #display message
  message("Plotting Weekly Dotplots...")
  
  #get names of metrics for plotting
  mets = colnames(full.metics)[-c(1:3)]
  
  for (i in seq_along(groups)) {
    
    #group for comparison and plotting
    group = groups[[i]]
    
    for (w in seq_along(weeks)) {
     
      #week for plotting
      week = weeks[w]
      
      for (j in seq_along(mets)) {
        
        #metric to plot
        met = mets[j]
        met.lab = gsub("\\.", " ", met)
        
        #filter data
        plot.df = filter(full.metics, Line %in% group & Week == week) %>% select(., Plate, Week, Line, met)
        colnames(plot.df)[4] <- "Metric"
        
        #filter stats
        if(length(group) == 2){
          plot.stats = filter(stats, Control.Group == group[1] & Comparison.Group == group[2] & Metric == met.lab & Week == week)
        } else if(length(group) == 3){
          plot.stats = filter(stats, Control.Group == group[1] & Comparison.Group %in% c(group[2], group[3]) & Metric == met.lab & Week == week)
        } else if(length(group) == 4){
          plot.stats = filter(stats, Control.Group == group[1] & Comparison.Group %in% c(group[2], group[3], group[4]) & Metric == met.lab & Week == week)
        }
        
        #clean up data for plotting
        plot.df$Line = factor(plot.df$Line, levels = group)
        plot.df$Week = paste("Week", plot.df$Week)
        
        #clean up stats for plotting
        plot.stats$Line = plot.stats$Comparison.Group
        plot.stats$Week = paste("Week", plot.stats$Week)
        
        #clean up plot label
        met.lab = gsub("Avg", "", met.lab)
        met.lab = gsub("Number of", "", met.lab)
        met.lab = gsub("Within", "in", met.lab)
        met.lab = gsub("Cross Correlation", "CC", met.lab)
        met.lab = gsub("Full Width at Half Height of", "Width of", met.lab)
        
        #plot
        p <- ggplot(plot.df, aes(y=Metric, x = Line, fill = Line)) + 
          geom_point(alpha = 0.7, size = 5, color = "black", shape = 21, position = position_jitter(width = 0.2)) + 
          geom_text(data = plot.stats, aes(x = Line, y = y.Offset, label = Significance.BH), size = 10) +
          stat_summary(fun = mean, fun.min = mean, fun.max = mean, geom = "crossbar", color = "black", width = 0.5, lwd = 0.7) +
          sfm +
          scm +
          #scale_y_continuous(trans='log10') +
          ylab(met.lab)+
          theme(axis.title.y = element_text(size = 20, color = "black", face = "bold"),
                axis.title.x =  element_text(size = 20, color = "black", face = "bold"),
                axis.text.y = element_text(size = 18, color = "black", face = "bold"),
                axis.text.x = element_text(size = 18, color = "black", face = "bold"),
                axis.ticks.x = element_blank(),
                panel.background = element_blank(),
                axis.line = element_line(color = "black", size =1),
                legend.position = "none",
          )
        
        #normalize y axis so it matches density plot
        p <- p + ylim(ggplot_build(p)$layout$panel_params[[1]]$y.range)
        
        
        #plot density plot
        p2 <- ggplot(plot.df, aes(y=Metric,fill = Line)) + 
          geom_density(stat = "density", orientation = "y", na.rm = TRUE, alpha = 0.5, color = "black", lwd = 0.5) + 
          sfm+
          ylim(ggplot_build(p)$layout$panel_params[[1]]$y.range)+
          theme(axis.title.y = element_blank(),
                axis.title.x =  element_text(size = 20, color = "white", face = "bold"),
                axis.text.y = element_blank(),
                axis.text.x = element_text(size = 20, color = "white", face = "bold"),
                axis.ticks.x = element_blank(),
                axis.ticks.y = element_blank(),
                axis.line.x = element_blank(),
                axis.line.y = element_blank(),
                panel.background = element_blank(),
                axis.line = element_line(color = "black", size =1),
                legend.position = "none",
          )
        
        #create directory and save
        path.name = "Plotting Output/Dotplots - Weekly/"
        dir.create(path.name, showWarnings = FALSE)
        path.name = paste(path.name, names(groups[i]), sep = "")
        dir.create(path.name, showWarnings = FALSE)
        path.name = paste(path.name, "/Week ", week, sep = "")
        dir.create(path.name, showWarnings = FALSE)
        path.name2 = paste(path.name, "/PNG ", week, sep = "")
        dir.create(path.name2, showWarnings = FALSE)
        
        ggsave(paste(path.name, "/", met.lab, ".pdf", sep = ""), plot = p, width = width, height = height)
        ggsave(paste(path.name2, "/", met.lab, ".png", sep = ""), plot = p, width = width, height = height)
        ggsave(paste(path.name, "/", met.lab, " density.pdf", sep = ""), plot = p2, width = 0.64, height = height)
        
      }
    }
  }
}

plot_weekly_ecdf <- function(full.metics, groups, weeks, width, height){
  
  #display message
  message("Plotting Weekly ecdf plots...")
  
  #get names of metrics for plotting
  mets = colnames(full.metics)[-c(1:3)]
  
  for (i in seq_along(groups)) {
    
    #group for comparison and plotting
    group = groups[[i]]
    
    for (w in seq_along(weeks)) {
      
      #week for plotting
      week = weeks[w]
      
      for (j in seq_along(mets)) {
        
        #metric to plot
        met = mets[j]
        met.lab = gsub("\\.", " ", met)
        
        #filter data
        plot.df = filter(full.metics, Line %in% group & Week == week) %>% select(., Plate, Week, Line, met)
        colnames(plot.df)[4] <- "Metric"
          
        #clean up data for plotting
        plot.df$Line = factor(plot.df$Line, levels = group)
        plot.df$Week = paste("Week", plot.df$Week)
        
        #clean up plot label
        met.lab = gsub("Avg", "", met.lab)
        met.lab = gsub("Number of", "", met.lab)
        met.lab = gsub("Within", "in", met.lab)
        met.lab = gsub("Cross Correlation", "CC", met.lab)
        met.lab = gsub("Full Width at Half Height of", "Width of", met.lab)
        
        #plot
        p <- ggplot(plot.df, aes(x=Metric, fill = Line, color = Line)) + 
          stat_ecdf(geom = "line", lwd = 1) +
          stat_ecdf(geom="point", alpha=0.5, size = 2.7, aes(fill = Line, color = Line), shape = 21)+
          sfm +
          scm +
          #scale_x_continuous(trans='log10') +
          ylab("Cumulative Probability") +
          xlab(met.lab)+
          theme(axis.title.y = element_text(size = 20, color = "black", face = "bold"),
                axis.title.x =  element_text(size = 20, color = "black", face = "bold"),
                axis.text.y = element_text(size = 18, color = "black", face = "bold"),
                axis.text.x = element_text(size = 18, color = "black", face = "bold"),
                panel.background = element_blank(),
                axis.line = element_line(color = "black", size =1),
                legend.position = "none",
          )
          
          #create directory and save
          path.name = "Plotting Output/ECDF Plots - Weekly/"
          dir.create(path.name, showWarnings = FALSE)
          path.name = paste(path.name, names(groups[i]), sep = "")
          dir.create(path.name, showWarnings = FALSE)
          path.name = paste(path.name, "/Week ", week, sep = "")
          dir.create(path.name, showWarnings = FALSE)
          path.name2 = paste(path.name, "/PNG ", week, sep = "")
          dir.create(path.name2, showWarnings = FALSE)
          
          ggsave(paste(path.name, "/", met.lab, ".pdf", sep = ""), plot = p, width = width, height = height)
          ggsave(paste(path.name2, "/", met.lab, ".png", sep = ""), plot = p, width = width, height = height)

        }
      }
    }
  }

plot_weekly_violins <- function(full.metics, wilcox.stats, groups, weeks){
  
  #display message
  message("Plotting Weekly Violin Plots...")
  
  #get names of metrics for plotting
  mets = colnames(full.metics)[-c(1:3)]
  
  for (i in seq_along(groups)) {
    
    #group for comparison and plotting
    group = groups[[i]]
    
    for (w in seq_along(weeks)) {
      
      #week for plotting
      week = weeks[w]
      
      for (j in seq_along(mets)) {
        
        #metric to plot
        met = mets[j]
        met.lab = gsub("\\.", " ", met)
        
        #filter data
        #plot.df = filter(full.metics, Line %in% group & Week == week) %>% select(., Plate, Week, Line, met)
        plot.df = filter(full.metics, Line %in% group & Week %in% week) %>% select(., Plate, Week, Line, met)
        colnames(plot.df)[4] <- "Var"
        
        #filter stats
        if(length(group) == 2){
          plot.stats = filter(wilcox.stats, Control.Group == group[1] & Comparison.Group == group[2] & Metric == met.lab & Week == week)
        } else if(length(group) == 3){
          plot.stats = filter(wilcox.stats, Control.Group == group[1] & Comparison.Group %in% c(group[2], group[3]) & Metric == met.lab & Week == week)
        }
        
        #clean up data for plotting
        plot.df$Line = factor(plot.df$Line, levels = group)
        plot.df$Week = paste("Week", plot.df$Week)
        
        #clean up stats for plotting
        plot.stats$Line = plot.stats$Comparison.Group
        plot.stats$Week = paste("Week", plot.stats$Week)

        #clean up plot label
        met.lab = gsub("Avg", "", met.lab)
        met.lab = gsub("Number of", "", met.lab)
        met.lab = gsub("Within", "in", met.lab)
        met.lab = gsub("Cross Correlation", "CC", met.lab)
        met.lab = gsub("Full Width at Half Height of", "Width of", met.lab)
        
        #plot
        ggplot(plot.df, aes(y=Var, x = Line, fill = Line)) + 
          geom_violin(trim = FALSE, size = 1, scale = "width", color = 'black') + 
          #geom_violin(trim = FALSE, size = 1, color = 'black') + 
          geom_text(data = plot.stats, aes(x = Line, y = y.Offset, label = Significance), size = 10) +
          stat_summary(fun = mean, fun.min = mean, fun.max = mean, geom = "crossbar", color = "black", width = 0.5, lwd = 0.7) +
          sfm +
          scm +
          ylab(met.lab)+
          theme(axis.title.y = element_blank(),
                #axis.title.y = element_text(size = 20, color = "black", face = "bold")
                axis.title.x =  element_text(size = 20, color = "white", face = "bold"),
                axis.text.y = element_text(size = 20, color = "black", face = "bold"),
                axis.text.x = element_text(size = 20, color = "black", face = "bold"),
                axis.ticks.x = element_blank(),
                panel.background = element_blank(),
                axis.line = element_line(color = "black", size =1),
                legend.position = "none",
          )
        
        #create directory and save
        path.name = "Plotting Output/Violin Plots - Weekly/"
        dir.create(path.name, showWarnings = FALSE)
        path.name = paste(path.name, names(groups[i]), sep = "")
        dir.create(path.name, showWarnings = FALSE)
        path.name = paste(path.name, "/Week ", week, sep = "")
        dir.create(path.name, showWarnings = FALSE)
        
        ggsave(paste(path.name, "/", met.lab, ".pdf", sep = ""), width = 3.5, height = 5)
        
      }
    }
  }
  
}


plot_pooled_lines <- function(metric.averages, groups, weeks, stats, error.type, error.viz, width, height){
  
  #display message
  message("Plotting Pooled Line Plots...")
  
  #subset weeks to plot
  tidy.averages = filter(metric.averages, Week %in% weeks)
  
  for (i in seq_along(groups)) {
    
    #subset group for comparison and plotting
    group = groups[[i]]
    group.df = filter(tidy.averages, Line %in% group)
    
    #get metrics for plotting
    metrics = unique(group.df$Metric)
    
    for (j in seq_along(metrics)) {
      
      #metric to plot
      met = metrics[j]
      met.lab = gsub("\\.", " ", met)
      
      #filter data
      plot.df = filter(group.df, Metric == met)
      
      #filter stats
      if(length(group) == 2){
        plot.stats = filter(stats, Control.Group == group[1] & Comparison.Group == group[2] & Metric == met.lab & Week %in% weeks)
      } else if(length(group) == 3){
        plot.stats = filter(stats, Control.Group == group[1] & Comparison.Group %in% c(group[2], group[3]) & Metric == met.lab & Week %in% weeks)
      } else if(length(group) == 4){
        plot.stats = filter(stats, Control.Group == group[1] & Comparison.Group %in% c(group[2], group[3], group[4]) & Metric == met.lab & Week %in% weeks)
      }
      
      #clean up plot label
      met.lab = gsub("Avg", "", met.lab)
      met.lab = gsub("Number of", "", met.lab)
      met.lab = gsub("Within", "in", met.lab)
      met.lab = gsub("Cross Correlation", "CC", met.lab)
      met.lab = gsub("Full Width at Half Height of", "Width of", met.lab)
      
      #clean up data for plotting
      plot.df$Line = factor(plot.df$Line, levels = group)
      plot.df$Week = as.numeric(plot.df$Week)
      
      #clean up stats for plotting
      plot.stats$Line = plot.stats$Comparison.Group
      plot.stats$Week = as.numeric(plot.stats$Week)
      plot.stats$y.Offset = by(plot.df$SE.ypos, list(plot.df$Week), max) + mean(plot.df$Value.Mean)*0.1
      
      
      
      ##############################################
      #     plotting with confidence intervals     #
      ##############################################
      if(error.type == "CI"){
        
        
        #Visualize with error bars
        if(error.viz == "errorbar"){
          #plot
          ggplot(plot.df, aes(y=Value.Mean, x = Week, group = Line)) + 
            geom_line(aes(color = Line), lwd = 1) +
            geom_point(aes(color = Line), size = 2)+
            geom_errorbar(aes(ymin=CI.yneg, ymax=CI.ypos, color = Line), lwd=1, width=0.15)+
            geom_text(data = plot.stats, aes(x = Week, y = y.Offset, label = Significance.BH), size = 10) +
            sfm +
            scm +
            ylab(met.lab)+
            xlab("Week")+
            xlim(min(weeks), max(weeks)) +
            theme(axis.title.y = element_text(size = 18, color = "black", face = "bold"),
                  axis.title.x = element_text(size = 18, color = "black", face = "bold"),
                  axis.text.y = element_text(size = 16, color = "black", face = "bold"),
                  axis.text.x = element_text(size = 16, color = "black", face = "bold"),
                  panel.background = element_blank(),
                  axis.line = element_line(color = "black", size =1),
                  axis.ticks = element_line(color = "black", size =1),
                  legend.position = "none",
                  panel.spacing.x = unit(2, "lines")
            )
          #create directory and save
          path.name = "Plotting Output/Line Plots - Pooled CI/"
          dir.create(path.name, showWarnings = FALSE)
          path.name = paste(path.name, names(groups[i]), sep = "")
          dir.create(path.name, showWarnings = FALSE)
          path.name2 = paste(path.name, "/PNG ", sep = "")
          dir.create(path.name2, showWarnings = FALSE)
          ggsave(paste(path.name, "/", met.lab, ".pdf", sep = ""), width = width, height = height)
          ggsave(paste(path.name2, "/", met.lab, ".png", sep = ""), width = width, height = height)
          
        }
        
        
        #Visualize with shaded region
        if(error.viz == "shade"){
          #plot
          ggplot(plot.df, aes(y=Value.Mean, x = Week, group = Line)) + 
            geom_ribbon(aes(ymin = CI.yneg, ymax = CI.ypos, fill = Line), alpha = 0.4) + 
            geom_line(aes(color = Line), lwd = 1) +
            geom_point(aes(color = Line), size = 2)+
            geom_text(data = plot.stats, aes(x = Week, y = y.Offset, label = Significance.BH), size = 10) +
            sfm +
            scm +
            ylab(met.lab)+
            xlab("Week")+
            xlim(min(weeks), max(weeks)) +
            theme(axis.title.y = element_text(size = 18, color = "black", face = "bold"),
                  axis.title.x = element_text(size = 18, color = "black", face = "bold"),
                  axis.text.y = element_text(size = 16, color = "black", face = "bold"),
                  axis.text.x = element_text(size = 16, color = "black", face = "bold"),
                  panel.background = element_blank(),
                  axis.line = element_line(color = "black", size =1),
                  axis.ticks = element_line(color = "black", size =1),
                  legend.position = "none",
                  panel.spacing.x = unit(2, "lines")
            )
          #create directory and save
          path.name = "Plotting Output/Line Plots - Pooled CI/"
          dir.create(path.name, showWarnings = FALSE)
          path.name = paste(path.name, names(groups[i]), sep = "")
          dir.create(path.name, showWarnings = FALSE)
          path.name2 = paste(path.name, "/PNG ", sep = "")
          dir.create(path.name2, showWarnings = FALSE)
          ggsave(paste(path.name, "/", met.lab, ".pdf", sep = ""), width = width, height = height)
          ggsave(paste(path.name2, "/", met.lab, ".png", sep = ""), width = width, height = height)
          
        }
      }
      
      
      ##############################################
      #        plotting with standard error        #
      ##############################################
      if(error.type == "SEM"){
        
        
        #Visualize with error bars
        if(error.viz == "errorbar"){
          #plot
          ggplot(plot.df, aes(y=Value.Mean, x = Week, group = Line)) + 
            geom_line(aes(color = Line), lwd = 1) +
            geom_point(aes(color = Line), size = 2)+
            geom_errorbar(aes(ymin=SE.yneg, ymax=SE.ypos, color = Line), lwd=1, width=0.15)+
            geom_text(data = plot.stats, aes(x = Week, y = y.Offset, label = Significance.BH), size = 10) +
            sfm +
            scm +
            ylab(met.lab)+
            xlab("Week")+
            xlim(min(weeks), max(weeks)) +
            theme(axis.title.y = element_text(size = 18, color = "black", face = "bold"),
                  axis.title.x = element_text(size = 18, color = "black", face = "bold"),
                  axis.text.y = element_text(size = 16, color = "black", face = "bold"),
                  axis.text.x = element_text(size = 16, color = "black", face = "bold"),
                  panel.background = element_blank(),
                  axis.line = element_line(color = "black", size =1),
                  axis.ticks = element_line(color = "black", size =1),
                  legend.position = "none",
                  panel.spacing.x = unit(2, "lines")
            )
          #create directory and save
          path.name = "Plotting Output/Line Plots - Pooled SEM/"
          dir.create(path.name, showWarnings = FALSE)
          path.name = paste(path.name, names(groups[i]), sep = "")
          dir.create(path.name, showWarnings = FALSE)
          path.name2 = paste(path.name, "/PNG ", sep = "")
          dir.create(path.name2, showWarnings = FALSE)
          ggsave(paste(path.name, "/", met.lab, ".pdf", sep = ""), width = width, height = height)
          ggsave(paste(path.name2, "/", met.lab, ".png", sep = ""), width = width, height = height)
        }
        
        
        #Visualize with shaded region
        if(error.viz == "shade"){
          ggplot(plot.df, aes(y=Value.Mean, x = Week, group = Line)) + 
            geom_ribbon(aes(ymin = SE.yneg, ymax = SE.ypos, fill = Line), alpha = 0.4) + 
            geom_line(aes(color = Line), lwd = 1) +
            geom_point(aes(color = Line), size = 2)+
            geom_text(data = plot.stats, aes(x = Week, y = y.Offset, label = Significance.BH), size = 10) +
            sfm +
            scm +
            ylab(met.lab)+
            xlab("Week")+
            xlim(min(weeks), max(weeks)) +
            theme(axis.title.y = element_text(size = 18, color = "black", face = "bold"),
                  axis.title.x = element_text(size = 18, color = "black", face = "bold"),
                  axis.text.y = element_text(size = 16, color = "black", face = "bold"),
                  axis.text.x = element_text(size = 16, color = "black", face = "bold"),
                  panel.background = element_blank(),
                  axis.line = element_line(color = "black", size =1),
                  axis.ticks = element_line(color = "black", size =1),
                  legend.position = "none",
                  panel.spacing.x = unit(2, "lines")
            )
          #create directory and save
          path.name = "Plotting Output/Line Plots - Pooled SEM/"
          dir.create(path.name, showWarnings = FALSE)
          path.name = paste(path.name, names(groups[i]), sep = "")
          dir.create(path.name, showWarnings = FALSE)
          path.name2 = paste(path.name, "/PNG ", sep = "")
          dir.create(path.name2, showWarnings = FALSE)
          ggsave(paste(path.name, "/", met.lab, ".pdf", sep = ""), width = width, height = height)
          ggsave(paste(path.name2, "/", met.lab, ".png", sep = ""), width = width, height = height)
        }
      }
    }
  }
}

plot_lines_facet <- function(metric.averages, weeks, stats, error.type, error.viz, width, height){
  
  #display message
  message("Plotting Pooled Line Plots...")
  
  #subset weeks to plot
  dat = filter(metric.averages, Week %in% weeks & Line %in% c("R841X", "R841X-C", "R841X DHPG", "R841X-C DHPG"))
  
  #add facet group
  dat$Group <- ""
  dat$Group[which(dat$Line %in% c("R841X", "R841X DHPG"))] <- "R841X"
  dat$Group[which(dat$Line %in% c("R841X-C", "R841X-C DHPG"))] <- "R841X-C"
  
  #get metrics for plotting
  metrics = unique(dat$Metric)
  
  for (j in seq_along(metrics)) {
    
    #metric to plot
    met = metrics[j]
    met.lab = gsub("\\.", " ", met)
    
    #filter data
    plot.df = filter(dat, Metric == met)
    
    #filter stats
    s1 = filter(stats, Control.Group == "R841X" & Comparison.Group == "R841X DHPG" & Metric == met.lab & Week %in% weeks)
    s1$Group = "R841X"
    s2 = filter(stats, Control.Group == "R841X-C" & Comparison.Group == "R841X-C DHPG" & Metric == met.lab & Week %in% weeks)
    s2$Group = "R841X-C"
    plot.stats = rbind(s1,s2)
    
    #clean up plot label
    met.lab = gsub("Avg", "", met.lab)
    met.lab = gsub("Number of", "", met.lab)
    met.lab = gsub("Within", "in", met.lab)
    met.lab = gsub("Cross Correlation", "CC", met.lab)
    met.lab = gsub("Full Width at Half Height of", "Width of", met.lab)
    
    #clean up data for plotting
    plot.df$Line = factor(plot.df$Line, levels = group)
    plot.df$Week = as.numeric(plot.df$Week)
    
    #clean up stats for plotting
    #plot.stats$Line = plot.stats$Comparison.Group
    plot.stats$Week = as.numeric(plot.stats$Week)
    plot.stats$y.Offset = by(plot.df$SE.ypos, list(plot.df$Week), max) + mean(plot.df$Value.Mean)*0.1
    
    
    
    ##############################################
    #     plotting with confidence intervals     #
    ##############################################
    if(error.type == "CI"){
      
      
      #Visualize with error bars
      if(error.viz == "errorbar"){
        #plot
        ggplot(plot.df, aes(y=Value.Mean, x = Week, group = Line)) + 
          geom_line(aes(color = Line), lwd = 1) +
          geom_point(aes(color = Line), size = 2)+
          geom_errorbar(aes(ymin=CI.yneg, ymax=CI.ypos, color = Line), lwd=1, width=0.15)+
          geom_text(data = plot.stats, aes(x = Week, y = y.Offset, label = Significance.BH), size = 10) +
          sfm +
          scm +
          ylab(met.lab)+
          xlab("Week")+
          xlim(min(weeks), max(weeks)) +
          theme(axis.title.y = element_text(size = 18, color = "black", face = "bold"),
                axis.title.x = element_text(size = 18, color = "black", face = "bold"),
                axis.text.y = element_text(size = 16, color = "black", face = "bold"),
                axis.text.x = element_text(size = 16, color = "black", face = "bold"),
                panel.background = element_blank(),
                axis.line = element_line(color = "black", size =1),
                axis.ticks = element_line(color = "black", size =1),
                legend.position = "none",
                panel.spacing.x = unit(2, "lines")
          )
        #create directory and save
        path.name = "Plotting Output/Line Plots - Pooled CI/"
        dir.create(path.name, showWarnings = FALSE)
        path.name = paste(path.name, names(groups[i]), sep = "")
        dir.create(path.name, showWarnings = FALSE)
        path.name2 = paste(path.name, "/PNG ", sep = "")
        dir.create(path.name2, showWarnings = FALSE)
        ggsave(paste(path.name, "/", met.lab, ".pdf", sep = ""), width = width, height = height)
        ggsave(paste(path.name2, "/", met.lab, ".png", sep = ""), width = width, height = height)
        
      }
      
      
      #Visualize with shaded region
      if(error.viz == "shade"){
        #plot
        ggplot(plot.df, aes(y=Value.Mean, x = Week, group = Line)) + 
          geom_ribbon(aes(ymin = CI.yneg, ymax = CI.ypos, fill = Line), alpha = 0.4) + 
          geom_line(aes(color = Line), lwd = 1) +
          geom_point(aes(color = Line), size = 2)+
          geom_text(data = plot.stats, aes(x = Week, y = y.Offset, label = Significance.BH), size = 10) +
          sfm +
          scm +
          ylab(met.lab)+
          xlab("Week")+
          xlim(min(weeks), max(weeks)) +
          theme(axis.title.y = element_text(size = 18, color = "black", face = "bold"),
                axis.title.x = element_text(size = 18, color = "black", face = "bold"),
                axis.text.y = element_text(size = 16, color = "black", face = "bold"),
                axis.text.x = element_text(size = 16, color = "black", face = "bold"),
                panel.background = element_blank(),
                axis.line = element_line(color = "black", size =1),
                axis.ticks = element_line(color = "black", size =1),
                legend.position = "none",
                panel.spacing.x = unit(2, "lines")
          )
        #create directory and save
        path.name = "Plotting Output/Line Plots - Pooled CI/"
        dir.create(path.name, showWarnings = FALSE)
        path.name = paste(path.name, names(groups[i]), sep = "")
        dir.create(path.name, showWarnings = FALSE)
        path.name2 = paste(path.name, "/PNG ", sep = "")
        dir.create(path.name2, showWarnings = FALSE)
        ggsave(paste(path.name, "/", met.lab, ".pdf", sep = ""), width = width, height = height)
        ggsave(paste(path.name2, "/", met.lab, ".png", sep = ""), width = width, height = height)
        
      }
    }
    
    
    ##############################################
    #        plotting with standard error        #
    ##############################################
    if(error.type == "SEM"){
      
      
      #Visualize with error bars
      if(error.viz == "errorbar"){
        #plot
        ggplot(plot.df, aes(y=Value.Mean, x = Week, group = Line)) + 
          geom_line(aes(color = Line), lwd = 1) +
          geom_point(aes(color = Line), size = 2)+
          geom_errorbar(aes(ymin=SE.yneg, ymax=SE.ypos, color = Line), lwd=1, width=0.15)+
          geom_text(data = plot.stats, aes(x = Week, y = y.Offset, label = Significance.BH), size = 10) +
          sfm +
          scm +
          ylab(met.lab)+
          xlab("Week")+
          xlim(min(weeks), max(weeks)) +
          theme(axis.title.y = element_text(size = 18, color = "black", face = "bold"),
                axis.title.x = element_text(size = 18, color = "black", face = "bold"),
                axis.text.y = element_text(size = 16, color = "black", face = "bold"),
                axis.text.x = element_text(size = 16, color = "black", face = "bold"),
                panel.background = element_blank(),
                axis.line = element_line(color = "black", size =1),
                axis.ticks = element_line(color = "black", size =1),
                legend.position = "none",
                panel.spacing.x = unit(2, "lines")
          )
        #create directory and save
        path.name = "Plotting Output/Line Plots - Pooled SEM/"
        dir.create(path.name, showWarnings = FALSE)
        path.name = paste(path.name, names(groups[i]), sep = "")
        dir.create(path.name, showWarnings = FALSE)
        path.name2 = paste(path.name, "/PNG ", sep = "")
        dir.create(path.name2, showWarnings = FALSE)
        ggsave(paste(path.name, "/", met.lab, ".pdf", sep = ""), width = width, height = height)
        ggsave(paste(path.name2, "/", met.lab, ".png", sep = ""), width = width, height = height)
      }
      
      
      #Visualize with shaded region
      if(error.viz == "shade"){
        ggplot(plot.df, aes(y=Value.Mean, x = Week, group = Line)) + 
          facet_wrap(.~ Group) +
          geom_ribbon(aes(ymin = SE.yneg, ymax = SE.ypos, fill = Line), alpha = 0.4) + 
          geom_line(aes(color = Line), lwd = 1) +
          geom_point(aes(color = Line), size = 2)+
          geom_text(data = plot.stats, aes(x = Week, y = y.Offset, label = Significance.BH), size = 10) +
          sfm +
          scm +
          ylab(met.lab)+
          xlab("Week")+
          xlim(min(weeks), max(weeks)) +
          theme(axis.title.y = element_text(size = 18, color = "black", face = "bold"),
                axis.title.x = element_text(size = 18, color = "black", face = "bold"),
                axis.text.y = element_text(size = 16, color = "black", face = "bold"),
                axis.text.x = element_text(size = 16, color = "black", face = "bold"),
                panel.background = element_blank(),
                axis.line = element_line(color = "black", size =1),
                axis.ticks = element_line(color = "black", size =1),
                legend.position = "none",
                panel.spacing.x = unit(2, "lines")
          )
        #create directory and save
        path.name = "Plotting Output/Line Plots - Pooled SEM/"
        dir.create(path.name, showWarnings = FALSE)
        path.name = paste(path.name, names(groups[i]), sep = "")
        dir.create(path.name, showWarnings = FALSE)
        path.name2 = paste(path.name, "/PNG ", sep = "")
        dir.create(path.name2, showWarnings = FALSE)
        ggsave(paste(path.name, "/", met.lab, ".pdf", sep = ""), width = width, height = height)
        ggsave(paste(path.name2, "/", met.lab, ".png", sep = ""), width = width, height = height)
      }
    }
  }
}




##### CorSE_Functions #####
Compile_CorSE <- function(bd, path){
  
  path = "G:/Synchrony Analysis/CorSE Analysis 2023"
  corfiles = list.files(path = path, pattern = "CorrData.csv", recursive = TRUE)
  
  cor.df = data.frame()
  for (i in 1:length(corfiles)){
    
    #get metadata
    corfile = corfiles[i]
    fname = str_split(corfile, pattern = "/") %>% unlist(.) %>% .[3]
    plate = str_split(fname, pattern = "_") %>% unlist(.) %>% .[1] %>% gsub("Plate-", "", .) %>% as.numeric(.)
    day = str_split(fname, pattern = "_") %>% unlist(.) %>% .[2] %>% gsub("D", "", .) %>% as.numeric(.)
    well = str_split(fname, pattern = "_") %>% unlist(.) %>% .[3]
    line = str_split(fname, pattern = "_") %>% unlist(.) %>% .[4]
    week = day2week_v2(day)
    
    #metaadata
    temp.df = data.frame(Plate = plate,
                         Day = day,
                         Week = week,
                         Well = well,
                         Line = line)
    
    #load data
    cordata = read_csv(paste(path,corfile,sep="/"), col_names = FALSE, show_col_types = FALSE) %>% unlist(.) %>% unname(.)
    
    #mean
    meancor = mean(cordata[-c(which(cordata == 0))])
    temp.df$Mean_Cor = meancor
    
    #strong connections
    thresholds = seq(0.25,0.75, by = 0.05)
    for (t in 1:length(thresholds)) {
      
      thresh = thresholds[t]
      n.strong = length(which(cordata[-c(which(cordata == 0))] > thresh))
      temp.s = data.frame(S = n.strong)
      colnames(temp.s) <- paste("n_Strong", thresh, sep = "_")
      temp.df = cbind(temp.df, temp.s)
      
    }
    
    #bind cor values
    temp.df = cbind(temp.df,as.data.frame(t(cordata)))
    
    #append
    cor.df = rbind(cor.df, temp.df)
    
  }
  
  dir.create(paste(bd, "CorSE Analysis", sep = "/"), showWarnings = FALSE)
  write.csv(cor.df, paste(bd, "/CorSE Analysis/", "Compiled_CorSE_data.csv", sep = ""), row.names = FALSE)
  write.csv(cor.df[,c(1:17)], paste(bd, "/CorSE Analysis/", "Compiled_CorSE_stats.csv", sep = ""), row.names = FALSE)
  
  
}

CorSE_Map_Defined <- function(cordat, maxChannels, genotype, cor.filename){
  
  #extract additional metadata 
  plate = meta.df$Plate[1]
  genotype = meta.df$Line[1]
  well = meta.df$Well[1]
  week = meta.df$Week[1]
  
  #set electrode labels
  elecs = as.character(c('11','21','31','41','51','61','71','81',
                         '12','22','32','42','52','62','72','82',
                         '13','23','33','43','53','63','73','83',
                         '14','24','34','44','54','64','74','84',
                         '15','25','35','45','55','65','75','85',
                         '16','26','36','46','56','66','76','86',
                         '17','27','37','47','57','67','77','87',
                         '18','28','38','48','58','68','78','88'
  ))
  
  #set up channel map
  nodes = data.frame(id = elecs,
                     x = rep(seq(1,8,by=1), 8),
                     y = c(rep(8,8), rep(7,8), rep(6,8), rep(5,8), rep(4,8), rep(3,8), rep(2,8), rep(1,8))
  )
  
  #Sort correlation values in descending order
  sortedvals = sort(as.numeric(unlist(cordat)), decreasing = TRUE)
  
  #Remove duplicate values
  sortedvals = sortedvals[seq(1,length(sortedvals),by=2)]
  sortedvals = unique(sortedvals)
  
  #initialize vectors
  source = vector()
  target = vector()
  weight = vector()
  
  #find indicies of cells that match top N connections
  for (i in 1:maxChannels) {
    
    idx = which(cordat == sortedvals[i], arr.ind = TRUE)
    
    s = elecs[idx[1,1]]
    t = elecs[idx[1,2]]
    w = sortedvals[i]
    
    source = c(source, s)
    target= c(target, t)
    weight = c(weight, w)
    
  }
  
  
  #define connections
  links = data.frame(from = source,
                     to = target,
                     weight = weight)
  
  #create network object
  mea.net <- graph_from_data_frame(d=links, vertices=nodes, directed=F) 
  
  #make edge widths scale with CorSE scores
  E(mea.net)$width <- E(mea.net)$weight*10
  
  #get color for plotting
  col = switch(genotype,
               "R841X"        = "#7FC97F",
               "R841X-C"      = "#FDAE61",
               "R841X DHPG"   = "#D1FFD1",
               "R841X-C DHPG" = "#FFDEBD",
               "CTRL"         = "#386CB0",
               "KO"           = "#F0027F",
  )
  
  #subset nodes to include only those involved in connections
  sub.nodes = filter(nodes, id %in% c(source, target))
  connected.nodes = which(V(mea.net)$name %in% sub.nodes$id)
  
  #set node colours
  col.vals = rep(rgb(.8, .8, .8, alpha = 0.3),length(V(mea.net)$name))
  col.vals[connected.nodes] = col
  
  #set node labels
  lab.vals = rep(NA,length(V(mea.net)$name))
  lab.vals[connected.nodes] = sub.nodes$id
  
  #set node label size
  lab.size = rep(NA,length(V(mea.net)$name))
  lab.size[connected.nodes] = 3
  
  #set node label colour
  lab.col = rep(NA,length(V(mea.net)$name))
  lab.col[connected.nodes] = "black"
  
  #assign color & label settings set in code above
  V(mea.net)$color = col.vals
  V(mea.net)$frame.color = col.vals
  V(mea.net)$label = lab.vals
  V(mea.net)$label.cex = lab.size
  V(mea.net)$label.color = lab.col
  
  
  #save plot as png
  #map.filename = gsub("CorrData.csv","CorSE-Map-defined.png",cor.filename)
  #map.filename = paste("Defined Maps - PNG", map.filename, sep = "\\")
  
  #png(map.filename, width = 1750, height = 1750)
  
  #plot(mea.net,
  #     edge.color = "black",
  #     edge.curved = 0.15)
  
  #dev.off()
  
  
  #resizing node labels for pdf plot
  lab.size = rep(NA,length(V(mea.net)$name))
  lab.size[connected.nodes] = 1
  V(mea.net)$label.cex = lab.size
  
  #save plot as pdf
  map.filename = gsub("CorrData.csv","CorSE-Map-defined.pdf",cor.filename)
  map.filename = paste(cpath, map.filename, sep = "/")
  
  pdf(map.filename, width = 7, height = 7)
  
  plot(mea.net,
       edge.color = "black",
       edge.curved = 0.15)
  
  dev.off()
  
  #save plot as object to return
  par(bg = 'white', mar = c(0,0,0,0), pty = "s")
  plt <- base2grob(~plot(mea.net, edge.color = "black", edge.curved = 0.15))
  
  return.list = list(plt, mea.net)
  
  return(return.list)
  
}

CorSE_Map_Threshold <- function(cordat, maxChannels, genotype, cor.filename){
  
  #extract additional metadata from CorrData filename
  split.filename = unlist(strsplit(cor.filename, "_"))
  plate = split.filename[1]
  genotype = split.filename[2]
  well = split.filename[3]
  week = split.filename[4]
  
  #set electrode labels
  elecs = as.character(c('11','21','31','41','51','61','71','81',
                         '12','22','32','42','52','62','72','82',
                         '13','23','33','43','53','63','73','83',
                         '14','24','34','44','54','64','74','84',
                         '15','25','35','45','55','65','75','85',
                         '16','26','36','46','56','66','76','86',
                         '17','27','37','47','57','67','77','87',
                         '18','28','38','48','58','68','78','88'
  ))
  
  
  #set up channel map
  nodes = data.frame(id = elecs,
                     x = rep(seq(1,8,by=1), 8),
                     y = c(rep(8,8), rep(7,8), rep(6,8), rep(5,8), rep(4,8), rep(3,8), rep(2,8), rep(1,8))
  )
  
  #Sort correlation values in descending order
  sortedvals = sort(as.numeric(unlist(cordat)), decreasing = TRUE)
  
  #Remove duplicate values
  sortedvals = sortedvals[seq(1,length(sortedvals),by=2)]
  sortedvals = unique(sortedvals)
  
  #subset values meeting threshold cutoff
  strongvals = sortedvals[which(sortedvals >= threshold)]
  
  #initialize vectors
  source = vector()
  target = vector()
  weight = vector()
  
  #find indicies of cells that match top N connections
  for (i in seq_along(strongvals)) {
    
    idx = which(cordat == strongvals[i], arr.ind = TRUE)
    
    s = elecs[idx[1,1]]
    t = elecs[idx[1,2]]
    w = sortedvals[i]
    
    source = c(source, s)
    target= c(target, t)
    weight = c(weight, w)
    
  }
  
  #define connections
  links = data.frame(from = source,
                     to = target,
                     weight = weight)
  
  #create network object
  mea.net <- graph_from_data_frame(d=links, vertices=nodes, directed=F) 
  
  #make edge widths scale with CorSE scores
  E(mea.net)$width <- E(mea.net)$weight*10
  
  #get color for plotting
  col = switch(genotype,
               "R841X"        = "#7FC97F",
               "R841X-C"      = "#FDAE61",
               "R841X DHPG"   = "#D1FFD1",
               "R841X-C DHPG" = "#FFDEBD",
               "KO"         = "#386CB0",
               "CTRL"           = "#F0027F",
  )
  
  #subset nodes to include only those involved in connections
  sub.nodes = filter(nodes, id %in% c(source, target))
  connected.nodes = which(V(mea.net)$name %in% sub.nodes$id)
  
  #get number of connections passing threshold
  passing.thresh = length(sub.nodes)
  
  #set node colours
  col.vals = rep(rgb(.8, .8, .8, alpha = 0.3),length(V(mea.net)$name))
  col.vals[connected.nodes] = col
  
  #set node labels
  lab.vals = rep(NA,length(V(mea.net)$name))
  lab.vals[connected.nodes] = sub.nodes$id
  
  #set node label size
  lab.size = rep(NA,length(V(mea.net)$name))
  lab.size[connected.nodes] = 3
  
  #set node label colour
  lab.col = rep(NA,length(V(mea.net)$name))
  lab.col[connected.nodes] = "black"
  
  #assign color & label settings set in code above
  V(mea.net)$color = col.vals
  V(mea.net)$frame.color = col.vals
  V(mea.net)$label = lab.vals
  V(mea.net)$label.cex = lab.size
  V(mea.net)$label.color = lab.col
  
  
  #save plot as png
  #map.filename = gsub("CorrData.csv","CorSE-Map-Threshold.png",cor.filename)
  #map.filename = paste("Threshold Maps - PNG", map.filename, sep = "\\")
  
  #png(map.filename, width = 1750, height = 1750)
  
  #plot(mea.net,
  #     edge.color = "black",
  #     edge.curved = 0.15)
  
  #dev.off()
  
  
  #resizing node labels for pdf plot
  lab.size = rep(NA,length(V(mea.net)$name))
  lab.size[connected.nodes] = 1
  V(mea.net)$label.cex = lab.size
  
  #save plot as pdf
  #map.filename = gsub("CorrData.csv","CorSE-Map-Threshold.pdf",cor.filename)
  #map.filename = paste("Threshold Maps - PDF", map.filename, sep = "\\")
  
  pdf(map.filename, width = 7, height = 7)
  
  plot(mea.net,
       edge.color = "black",
       edge.curved = 0.15)
  
  dev.off()
  
  #save plot as object to return
  par(bg = 'white', mar = c(0,0,0,0), pty = "s")
  plt <- base2grob(~plot(mea.net,
                         edge.color = "black", 
                         edge.curved = 0.15))
  
  return.list = list(plt, mea.net, passing.thresh)
  
  return(return.list)
  
}

CorSE_Heatmap <- function(cordat, cor.filename, palette, direct, scale.min, scale.max, legend){
  
  #extract additional metadata from CorrData filename
  split.filename = unlist(strsplit(cor.filename, "_"))
  plate = split.filename[1]
  genotype = split.filename[2]
  well = split.filename[3]
  week = split.filename[4]
  
  #set electrode labels
  elecs = as.character(c('11','21','31','41','51','61','71','81',
                         '12','22','32','42','52','62','72','82',
                         '13','23','33','43','53','63','73','83',
                         '14','24','34','44','54','64','74','84',
                         '15','25','35','45','55','65','75','85',
                         '16','26','36','46','56','66','76','86',
                         '17','27','37','47','57','67','77','87',
                         '18','28','38','48','58','68','78','88'
  ))
  
  #convert cordat to matrix
  cormat = as.matrix(cordat)
  
  #add row and column names
  rownames(cormat) <- paste("t", elecs, sep="")
  colnames(cormat) <- paste("s", elecs, sep="")
  
  #reshape for ggplot
  meltdat = reshape2::melt(cormat)
  
  gradcols = cubehelix(100)
  
  #plot
  plt <- ggplot(data = meltdat, aes(x=Var1, y=Var2, fill=value)) + 
    geom_tile() +
    #scale_fill_scico("CorSE \n Score", palette = palette, direction = direct, midpoint = 0, limits = c(scale.min, scale.max)) +
    scale_fill_gradientn(colours = gradcols, limits=c(scale.min, scale.max), guide = "colourbar", na.value = gradcols[1]) +
    labs(title = paste(well, " (", genotype, ") - ", gsub("-", " ", week), sep = ""), x=NULL, y=NULL) + 
    theme(
      axis.ticks=element_blank(),
      axis.text=element_blank(),
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
      aspect.ratio = 1,
      legend.position = legend
    )
  
  plt
  
  #save as png
  heatmap.filename = gsub("CorrData.csv","CorSE_Heatmap.png",cor.filename)
  heatmap.filename = paste("Heatmaps - PNG", heatmap.filename, sep = "\\")
  ggsave(heatmap.filename, height = 7, width = 7)
  
  #save as png
  heatmap.filename = gsub("CorrData.csv","CorSE_Heatmap.pdf",cor.filename)
  heatmap.filename = paste("Heatmaps - PDF", heatmap.filename, sep = "\\")
  ggsave(heatmap.filename, height = 7, width = 7)
  
  #return variables
  return.list = list(plt, meltdat)
  return(return.list)
  
}

GIF.export <- function(plot.list, labels, filename, interval, width, height, hjust){
  
  if(missing(labels)){
    labels = rep("", length(plot.list))
  }
  
  if(missing(hjust)){
    hjust = -1
  }
  
  saveGIF({
    
    for (g in seq_along(plot.list)){
      
      lab = labels[g]
      patch = plot.list[[g]]
      newpatch <- plot_grid(patch, labels = lab, hjust = hjust)
      plot(newpatch)
      
    }
    
  },
  
  interval = interval,
  movie.name=filename,
  ani.width = width,
  ani.height = height,
  autoplay = FALSE)
  
}

fix.electrode.errors <- function(cordat, plate){
  
  #fixing error with plate 63 where electrode 16 is deactivated in some recordings
  if(plate == "Plate 43" & length(cordat) == 63){
    
    #padding with NA
    cordat$E64 = NA
    cordat[64,] = NA
    
    #getting indicies of deactivated electrode 16
    elec.idx = which(elecs == '16')
    idx.lead = elec.idx - 1
    
    #reshaping
    order.new = c(1:idx.lead, 64, elec.idx:63)
    cordat = cordat[order.new, order.new]
    
    #fixing diagonal
    cordat[elec.idx,elec.idx] = 0
    
    #fix column names
    colnames(cordat) <- paste("X", 1:64, sep = "")
    
  }
  
  return(cordat)
  
}

subset_weeks <- function(corfiles, weeks.include){
  
  sub.corfiles = vector()
  
  for (w in seq_along(weeks.include)) {
    
    w.files = corfiles[str_which(corfiles, weeks.include[w])]
    sub.corfiles = c(sub.corfiles, w.files)
    
  }
  
  return(sub.corfiles)
  
}

Sumarize_CorSE <- function(corfiles, pos.thresh, neg.thresh){
  
  CorSE.data = as.data.frame(matrix(ncol = 11, nrow = length(corfiles)))
  
  for (i in seq_along(corfiles)) {
    
    #load data
    cor.filename = corfiles[i]
    cordat <- read_csv(cor.filename, col_names = FALSE, show_col_types = FALSE)
    
    #get metadata
    genotype = unlist(strsplit(cor.filename, "/"))[4] %>% strsplit(., "_") %>% unlist() %>% .[2]
    plate = unlist(strsplit(cor.filename, "/"))[2] 
    plate.num = plate %>% gsub("Plate ", "", .) %>% as.numeric()
    well = unlist(strsplit(cor.filename, "/"))[3]
    week = unlist(strsplit(cor.filename, "/"))[4] %>% strsplit(., "_") %>% unlist() %>% .[4] %>% gsub("Week-", "", .) %>% as.numeric()
    
    #mean CorSE score
    mean.cor = mean(unlist(cordat))
    
    #mean CorSE score - positive correlations only
    mean.pos = unlist(cordat) %>% subset(unlist(cordat > 0)) %>% mean()
    
    #number of positive correlatons
    n.pos = unlist(cordat) %>% subset(unlist(cordat > 0)) %>% length()
    
    #number of stong positive correlatons
    n.pos.strong = unlist(cordat) %>% subset(unlist(cordat >= pos.thresh)) %>% length()
    
    #mean CorSE score - negative correlations only
    mean.neg = unlist(cordat) %>% subset(unlist(cordat < 0)) %>% mean()
    
    #number of negative correlatons
    n.neg = unlist(cordat) %>% subset(unlist(cordat < 0)) %>% length()
    
    #number of stong positive correlatons
    n.neg.strong = unlist(cordat) %>% subset(unlist(cordat <= neg.thresh)) %>% length()
    
    #creat summary data frame
    sum.df = data.frame(Plate = plate.num,
                        Week = week,
                        Well = well,
                        Genotype = genotype,
                        Mean.CorSE = mean.cor,
                        Mean.Pos = mean.pos,
                        Mean.Neg = mean.neg,
                        n.Pos.Strong = n.pos.strong,
                        n.Neg.Strong = n.neg.strong,
                        n.Pos.Total = n.pos,
                        n.Neg.Total = n.neg
    )
    
    #append results to master data frame
    CorSE.data[i,] = sum.df
    
  }
  
  #fix column names
  colnames(CorSE.data) <- colnames(sum.df)
  
  #save
  return(CorSE.data)
  
  
}

get_averages_CorSE <- function(dat){
  
  metric.averages = data.frame()
  lines = unique(dat$Line)
  mets = c("Mean_Cor", "n_Strong_0.5","n_Strong_0.55","n_Strong_0.6","n_Strong_0.65","n_Strong_0.7","n_Strong_0.75")
  
  
  for (i in seq_along(lines)) {
    
    #subset line
    line = lines[i]
    line.df = filter(dat, Line == line)
    weeks = unique(line.df$Week)
    
    for (j in seq_along(weeks)) {
      
      #subset week
      week = weeks[j]
      week.df = filter(line.df, Week == week)
      
      for (k in seq_along(mets)) {
        
        #subset metric
        metric = mets[k]
        met.df = select(week.df, all_of(metric))
        met.df = na.omit(met.df)
        colnames(met.df) <- "Metric"
        
        #skip if all NA
        if(nrow(met.df) == 0){
          next
        }
        
        #get mean and SE
        met.mean = mean(met.df$Metric)
        met.se = se(met.df$Metric)
        se.yneg = met.mean - met.se
        se.ypos = met.mean + met.se
        
        #get 95% confidence interval
        lmod <- lm(Metric ~ 1, met.df)
        ci <- confint(lmod, level=0.95)
        ci.yneg = ci[1]
        ci.ypos = ci[2]
        
        #create data frame for results
        temp.df = data.frame(Line = line,
                             Week = week,
                             Metric = metric,
                             Value.Mean = met.mean,
                             Value.SE = met.se,
                             SE.yneg = se.yneg,
                             SE.ypos = se.ypos,
                             CI.yneg = ci.yneg,
                             CI.ypos = ci.ypos)
        
        #append
        metric.averages = rbind(metric.averages, temp.df)
        
      }
    }
  }
  return(metric.averages)
}

plot_CorSE_lines <- function(avg.corse, groups, weeks, stats, error.type, error.viz, width, height){
  
  #display message
  message("Plotting Pooled Line Plots...")
  
  #subset weeks to plot
  avg.dat = filter(avg.corse, Week %in% weeks)
  
  for (i in seq_along(groups)) {
    
    #subset group for comparison and plotting
    group = groups[[i]]
    group.df = filter(avg.dat, Line %in% group)
    
    #get metrics for plotting
    metrics = unique(group.df$Metric)
    
    for (j in seq_along(metrics)) {
      
      #metric to plot
      met = metrics[j]
      met.lab = gsub("\\.", " ", met)
      
      #filter data
      plot.df = filter(group.df, Metric == met)
      
      #filter stats
      if(length(group) == 2){
        plot.stats = filter(stats, Control.Group == group[1] & Comparison.Group == group[2] & Metric == met.lab & Week %in% weeks)
      } else if(length(group) == 3){
        plot.stats = filter(stats, Control.Group == group[1] & Comparison.Group %in% c(group[2], group[3]) & Metric == met.lab & Week %in% weeks)
      } else if(length(group) == 4){
        plot.stats = filter(stats, Control.Group == group[1] & Comparison.Group %in% c(group[2], group[3], group[4]) & Metric == met.lab & Week %in% weeks)
      }
      
      #clean up data for plotting
      plot.df$Line = factor(plot.df$Line, levels = group)
      plot.df$Week = as.numeric(plot.df$Week)
      
      #clean up stats for plotting
      plot.stats$Line = plot.stats$Comparison.Group
      plot.stats$Week = as.numeric(plot.stats$Week)
      plot.stats$y.Offset = by(plot.df$CI.ypos, list(plot.df$Week), max) + mean(plot.df$Value.Mean)*0.1
      
      
      
      ##############################################
      #     plotting with confidence intervals     #
      ##############################################
      if(error.type == "CI"){
        
        
        #Visualize with error bars
        if(error.viz == "errorbar"){
          #plot
          ggplot(plot.df, aes(y=Value.Mean, x = Week, group = Line)) + 
            geom_line(aes(color = Line), lwd = 1) +
            geom_point(aes(color = Line), size = 2)+
            geom_errorbar(aes(ymin=CI.yneg, ymax=CI.ypos, color = Line), lwd=1, width=0.15)+
            geom_text(data = plot.stats, aes(x = Week, y = y.Offset, label = Significance), size = 10) +
            sfm +
            scm +
            ylab(met.lab)+
            xlab("Week")+
            xlim(min(weeks), max(weeks)) +
            theme(axis.title.y = element_text(size = 18, color = "black", face = "bold"),
                  axis.title.x = element_text(size = 18, color = "black", face = "bold"),
                  axis.text.y = element_text(size = 16, color = "black", face = "bold"),
                  axis.text.x = element_text(size = 16, color = "black", face = "bold"),
                  panel.background = element_blank(),
                  axis.line = element_line(color = "black", size =1),
                  axis.ticks = element_line(color = "black", size =1),
                  legend.position = "none",
                  panel.spacing.x = unit(2, "lines")
            )
          #create directory and save
          path.name = "C:/Users/Fraser/Desktop/SHANK2 Paper Data/Plotting Output/Line Plots - CorSE CI/"
          dir.create(path.name, showWarnings = FALSE)
          path.name = paste(path.name, names(groups[i]), sep = "")
          dir.create(path.name, showWarnings = FALSE)
          path.name2 = paste(path.name, "/PNG ", sep = "")
          dir.create(path.name2, showWarnings = FALSE)
          ggsave(paste(path.name, "/", met.lab, ".pdf", sep = ""), width = width, height = height)
          ggsave(paste(path.name2, "/", met.lab, ".png", sep = ""), width = width, height = height)
          
        }
        
        
        #Visualize with shaded region
        if(error.viz == "shade"){
          #plot
          ggplot(plot.df, aes(y=Value.Mean, x = Week, group = Line)) + 
            geom_ribbon(aes(ymin = CI.yneg, ymax = CI.ypos, fill = Line), alpha = 0.4) + 
            geom_line(aes(color = Line), lwd = 1) +
            geom_point(aes(color = Line), size = 2)+
            geom_text(data = plot.stats, aes(x = Week, y = y.Offset, label = Significance), size = 10) +
            sfm +
            scm +
            ylab(met.lab)+
            xlab("Week")+
            xlim(min(weeks), max(weeks)) +
            theme(axis.title.y = element_text(size = 18, color = "black", face = "bold"),
                  axis.title.x = element_text(size = 18, color = "black", face = "bold"),
                  axis.text.y = element_text(size = 16, color = "black", face = "bold"),
                  axis.text.x = element_text(size = 16, color = "black", face = "bold"),
                  panel.background = element_blank(),
                  axis.line = element_line(color = "black", size =1),
                  axis.ticks = element_line(color = "black", size =1),
                  legend.position = "none",
                  panel.spacing.x = unit(2, "lines")
            )
          #create directory and save
          path.name = "C:/Users/Fraser/Desktop/SHANK2 Paper Data/Plotting Output/Line Plots - CorSE CI/"
          dir.create(path.name, showWarnings = FALSE)
          path.name = paste(path.name, names(groups[i]), sep = "")
          dir.create(path.name, showWarnings = FALSE)
          path.name2 = paste(path.name, "/PNG ", sep = "")
          dir.create(path.name2, showWarnings = FALSE)
          ggsave(paste(path.name, "/", met.lab, ".pdf", sep = ""), width = width, height = height)
          ggsave(paste(path.name2, "/", met.lab, ".png", sep = ""), width = width, height = height)
          
        }
      }
      
      
      ##############################################
      #        plotting with standard error        #
      ##############################################
      if(error.type == "SEM"){
        
        
        #Visualize with error bars
        if(error.viz == "errorbar"){
          #plot
          ggplot(plot.df, aes(y=Value.Mean, x = Week, group = Line)) + 
            geom_line(aes(color = Line), lwd = 1) +
            geom_point(aes(color = Line), size = 2)+
            geom_errorbar(aes(ymin=SE.yneg, ymax=SE.ypos, color = Line), lwd=1, width=0.15)+
            geom_text(data = plot.stats, aes(x = Week, y = y.Offset, label = Significance), size = 10) +
            sfm +
            scm +
            ylab(met.lab)+
            xlab("Week")+
            xlim(min(weeks), max(weeks)) +
            theme(axis.title.y = element_text(size = 18, color = "black", face = "bold"),
                  axis.title.x = element_text(size = 18, color = "black", face = "bold"),
                  axis.text.y = element_text(size = 16, color = "black", face = "bold"),
                  axis.text.x = element_text(size = 16, color = "black", face = "bold"),
                  panel.background = element_blank(),
                  axis.line = element_line(color = "black", size =1),
                  axis.ticks = element_line(color = "black", size =1),
                  legend.position = "none",
                  panel.spacing.x = unit(2, "lines")
            )
          #create directory and save
          path.name = "C:/Users/Fraser/Desktop/SHANK2 Paper Data/Plotting Output/Line Plots - CorSE SEM/"
          dir.create(path.name, showWarnings = FALSE)
          path.name = paste(path.name, names(groups[i]), sep = "")
          dir.create(path.name, showWarnings = FALSE)
          path.name2 = paste(path.name, "/PNG ", sep = "")
          dir.create(path.name2, showWarnings = FALSE)
          ggsave(paste(path.name, "/", met.lab, ".pdf", sep = ""), width = width, height = height)
          ggsave(paste(path.name2, "/", met.lab, ".png", sep = ""), width = width, height = height)
        }
        
        
        #Visualize with shaded region
        if(error.viz == "shade"){
          ggplot(plot.df, aes(y=Value.Mean, x = Week, group = Line)) + 
            geom_ribbon(aes(ymin = SE.yneg, ymax = SE.ypos, fill = Line), alpha = 0.4) + 
            geom_line(aes(color = Line), lwd = 1) +
            geom_point(aes(color = Line), size = 2)+
            geom_text(data = plot.stats, aes(x = Week, y = y.Offset, label = Significance), size = 10) +
            sfm +
            scm +
            ylab(met.lab)+
            xlab("Week")+
            xlim(min(weeks), max(weeks)) +
            theme(axis.title.y = element_text(size = 18, color = "black", face = "bold"),
                  axis.title.x = element_text(size = 18, color = "black", face = "bold"),
                  axis.text.y = element_text(size = 16, color = "black", face = "bold"),
                  axis.text.x = element_text(size = 16, color = "black", face = "bold"),
                  panel.background = element_blank(),
                  axis.line = element_line(color = "black", size =1),
                  axis.ticks = element_line(color = "black", size =1),
                  legend.position = "none",
                  panel.spacing.x = unit(2, "lines")
            )
          #create directory and save
          path.name = "C:/Users/Fraser/Desktop/SHANK2 Paper Data/Plotting Output/Line Plots - CorSE SEM/"
          dir.create(path.name, showWarnings = FALSE)
          path.name = paste(path.name, names(groups[i]), sep = "")
          dir.create(path.name, showWarnings = FALSE)
          path.name2 = paste(path.name, "/PNG ", sep = "")
          dir.create(path.name2, showWarnings = FALSE)
          ggsave(paste(path.name, "/", met.lab, ".pdf", sep = ""), width = width, height = height)
          ggsave(paste(path.name2, "/", met.lab, ".png", sep = ""), width = width, height = height)
        }
      }
    }
  }
}

compile_stats_CorSE_AD <- function(dat, comparisons){
  
  mets = colnames(dat)[c(6:17)]
  sub.dat = select(dat, Week, Line, mets)
  AD.df = data.frame()
  
  #n for testing correction
  n = length(mets)
  
  for (i in seq_along(comparisons)) {
    
    lines = comparisons[[i]]
    ctrl.line = lines[1]
    comp.line = lines[2]
    
    group.df = filter(sub.dat, Line %in% lines)
    weeks = unique(group.df$Week)
    
    for (j in seq_along(weeks)) {
      
      #subset week
      week = weeks[j]
      week.df = filter(group.df, Week == week)
      
      AD.week = data.frame()
      
      for (k in seq_along(mets)){
        
        #subset metric
        met = mets[k]
        met.lab = gsub("\\.", " ", met)
        metric.df =  select(week.df, Week, Line, met)
        
        #subset control group data
        ctrl.dat = filter(metric.df, Line == ctrl.line)[[3]] %>% .[!is.na(.)]
        ctrl.n = length(ctrl.dat)
        ctrl.mean = mean(ctrl.dat)
        ctrl.se = sd(ctrl.dat) / sqrt(ctrl.n)
        
        #subset comparison group data
        comp.dat = filter(metric.df, Line == comp.line)[[3]] %>% .[!is.na(.)]
        comp.n = length(comp.dat)
        comp.mean = mean(comp.dat)
        comp.se = sd(comp.dat) / sqrt(comp.n)
        
        #Enter NS if no data
        if (ctrl.n == 0 | comp.n == 0){
          
          temp.df = data.frame(Week = week,
                               Metric = met.lab,
                               Control.Group = ctrl.line,
                               Comparison.Group = comp.line,
                               Control.n = ctrl.n,
                               Comparison.n = comp.n,
                               Control.Mean = NA,
                               Control.SEM = NA,
                               Comparison.Mean = NA,
                               Comparison.SEM = NA,
                               Fold.Change = NA,
                               y.Offset = NA,
                               AD.Statistic = NA,
                               p.value = NA,
                               Significance = NA,
                               p.adj.bon = NA,
                               Significance.Bon = NA)
          
        } else {
          
          #fold change and AD stats
          fc = comp.mean / ctrl.mean
          ad.res = ad.test(ctrl.dat, comp.dat)
          ad.stat = ad.res$ad[2,1]
          pval = ad.res$ad[2,3]
          
          #get significance star
          if (pval < 0.005){
            sig = "***"
          } else if (pval < 0.01){
            sig = "**"
          } else if (pval < 0.05){
            sig = "*"
          } else {
            sig = ""
          }
          
          
          #bonferroni correction
          p.adj.bon = p.adjust(pval, method = "bonferroni", n)
          
          #get significance star
          if (p.adj.bon < 0.005){
            sig.bon = "***"
          } else if (p.adj.bon < 0.01){
            sig.bon = "**"
          } else if (p.adj.bon < 0.05){
            sig.bon = "*"
          } else {
            sig.bon = ""
          }
          
          
          #y.offset = max(c(ctrl.dat, comp.dat)) + max(c(ctrl.dat, comp.dat)) * 0.2
          y.offset = max(c(mean(ctrl.dat), mean(comp.dat))) * 3
          
          
          #Create data frame
          temp.df = data.frame(Week = week,
                               Metric = met.lab,
                               Control.Group = ctrl.line,
                               Comparison.Group = comp.line,
                               Control.n = ctrl.n,
                               Comparison.n = comp.n,
                               Control.Mean = ctrl.mean,
                               Control.SEM = ctrl.se,
                               Comparison.Mean = comp.mean,
                               Comparison.SEM = comp.se,
                               Fold.Change = fc,
                               y.Offset = y.offset,
                               AD.Statistic = ad.stat,
                               p.value = pval,
                               Significance = sig,
                               p.adj.bon = p.adj.bon,
                               Significance.Bon = sig.bon
          )
          
        }
        #append
        AD.week = rbind(AD.week, temp.df)
      }
      
      #Apply BH correction by week
      AD.week$p.adj.BH = p.adjust(AD.week$p.value, method = "BH")
      AD.week$Significance.BH <- ""
      AD.week$Significance.BH[which(AD.week$p.adj.BH < 0.05)] <- "*"
      AD.week$Significance.BH[which(AD.week$p.adj.BH < 0.01)] <- "**"
      AD.week$Significance.BH[which(AD.week$p.adj.BH < 0.005)] <- "***"
      
      #append
      AD.df = rbind(AD.df, AD.week)
    }
  }
  write.csv(AD.df, "CorSE_AD_stats.csv")
  return(AD.df)
}

compile_stats_CorSE_WC <- function(dat, comparisons){
  
  mets = colnames(dat)[c(6:17)]
  sub.dat = select(dat, Week, Line, mets)
  AD.df = data.frame()
  
  #n for testing correction
  n = length(mets)
  
  for (i in seq_along(comparisons)) {
    
    lines = comparisons[[i]]
    ctrl.line = lines[1]
    comp.line = lines[2]
    
    group.df = filter(sub.dat, Line %in% lines)
    weeks = unique(group.df$Week)
    
    for (j in seq_along(weeks)) {
      
      #subset week
      week = weeks[j]
      week.df = filter(group.df, Week == week)
      
      AD.week = data.frame()
      
      for (k in seq_along(mets)){
        
        #subset metric
        met = mets[k]
        met.lab = gsub("\\.", " ", met)
        metric.df =  select(week.df, Week, Line, met)
        
        #subset control group data
        ctrl.dat = filter(metric.df, Line == ctrl.line)[[3]] %>% .[!is.na(.)]
        ctrl.n = length(ctrl.dat)
        ctrl.mean = mean(ctrl.dat)
        ctrl.se = sd(ctrl.dat) / sqrt(ctrl.n)
        
        #subset comparison group data
        comp.dat = filter(metric.df, Line == comp.line)[[3]] %>% .[!is.na(.)]
        comp.n = length(comp.dat)
        comp.mean = mean(comp.dat)
        comp.se = sd(comp.dat) / sqrt(comp.n)
        
        #Enter NS if no data
        if (ctrl.n == 0 | comp.n == 0){
          
          temp.df = data.frame(Week = week,
                               Metric = met.lab,
                               Control.Group = ctrl.line,
                               Comparison.Group = comp.line,
                               Control.n = ctrl.n,
                               Comparison.n = comp.n,
                               Control.Mean = NA,
                               Control.SEM = NA,
                               Comparison.Mean = NA,
                               Comparison.SEM = NA,
                               Fold.Change = NA,
                               y.Offset = NA,
                               AD.Statistic = NA,
                               p.value = NA,
                               Significance = NA,
                               p.adj.bon = NA,
                               Significance.Bon = NA)
          
        } else {
          
          #fold change and AD stats
          fc = comp.mean / ctrl.mean
          # ad.res = wilcox.test(ctrl.dat, comp.dat)
          # ad.stat = ad.res$ad[2,1]
          # pval = ad.res$ad[2,3]
          ad.res = wilcox.test(ctrl.dat, comp.dat)
          ad.stat = ad.res$statistic
          pval = ad.res$p.value
          
          #get significance star
          if (pval < 0.005){
            sig = "***"
          } else if (pval < 0.01){
            sig = "**"
          } else if (pval < 0.05){
            sig = "*"
          } else {
            sig = ""
          }
          
          
          #bonferroni correction
          p.adj.bon = p.adjust(pval, method = "bonferroni", n)
          
          #get significance star
          if (p.adj.bon < 0.005){
            sig.bon = "***"
          } else if (p.adj.bon < 0.01){
            sig.bon = "**"
          } else if (p.adj.bon < 0.05){
            sig.bon = "*"
          } else {
            sig.bon = ""
          }
          
          
          #y.offset = max(c(ctrl.dat, comp.dat)) + max(c(ctrl.dat, comp.dat)) * 0.2
          y.offset = max(c(mean(ctrl.dat), mean(comp.dat))) * 3
          
          
          #Create data frame
          temp.df = data.frame(Week = week,
                               Metric = met.lab,
                               Control.Group = ctrl.line,
                               Comparison.Group = comp.line,
                               Control.n = ctrl.n,
                               Comparison.n = comp.n,
                               Control.Mean = ctrl.mean,
                               Control.SEM = ctrl.se,
                               Comparison.Mean = comp.mean,
                               Comparison.SEM = comp.se,
                               Fold.Change = fc,
                               y.Offset = y.offset,
                               AD.Statistic = ad.stat,
                               p.value = pval,
                               Significance = sig,
                               p.adj.bon = p.adj.bon,
                               Significance.Bon = sig.bon
          )
          
        }
        #append
        AD.week = rbind(AD.week, temp.df)
      }
      
      #Apply BH correction by week
      AD.week$p.adj.BH = p.adjust(AD.week$p.value, method = "BH")
      AD.week$Significance.BH <- ""
      AD.week$Significance.BH[which(AD.week$p.adj.BH < 0.05)] <- "*"
      AD.week$Significance.BH[which(AD.week$p.adj.BH < 0.01)] <- "**"
      AD.week$Significance.BH[which(AD.week$p.adj.BH < 0.005)] <- "***"
      
      #append
      AD.df = rbind(AD.df, AD.week)
    }
  }
  write.csv(AD.df, "CorSE_AD_stats.csv")
  return(AD.df)
}


