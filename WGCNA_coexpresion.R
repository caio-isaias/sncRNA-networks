### LOAD PACKAGES
library(WGCNA)
library(edgeR)
library(variancePartition)
library(lme4)

### LOAD DATA

# processed data
miRNA_DGE = readRDS("data/DGEList_miRNA_filtered_norm.RDS")
piRNA_DGE = readRDS("data/DGEList_piRNA_filtered_norm.RDS")
tRNA_DGE  = readRDS("data/DGEList_tRNA_filtered_norm.RDS")
# sncRNA dynamic classification results
dynamic_candidates  = readRDS("data/dynamic_cadidate_all_subjects_DC_DE.RDS")

info = read.csv("data/data_for_modeling.csv",
                row.names = read.csv("data/data_for_modeling.csv")[,1])[,-1]

# filter out non-dynamic sncRNA
dynamic_miRNA = dynamic_candidates$miRNA_DE[(dynamic_candidates$miRNA_DE$dynamic_count_DE >= 7),]$miRNA
dynamic_piRNA = dynamic_candidates$piRNA_DE[(dynamic_candidates$piRNA_DE$dynamic_count_DE >= 7),]$piRNA
dynamic_tRNA = dynamic_candidates$tRNA_DE[(dynamic_candidates$tRNA_DE$dynamic_count_DE >= 7),]$tRNA

sncRNA_DGE = rbind(miRNA_DGE[dynamic_miRNA,rownames(info)],
                   piRNA_DGE[dynamic_piRNA,rownames(info)],
                   tRNA_DGE[dynamic_tRNA,rownames(info)])

# get the dynamic PSS reference visits
increase_dynamic = read.csv("data/mark_increase_dynamic.csv")[,-1]
stability_low_dynamic = read.csv("data/mark_stability_low_dynamic.csv")[,-1]
stability_moderate_dynamic = read.csv("data/mark_stability_moderate_dynamic.csv")[,-1]

# function to get samples from PSS dynamic classification data.frame  
get_samples <- function(mark_df, delay = 0){
    subjects = gsub("50","",mark_df$Record_ID) # remove "50", from the Record_ID
    visit = mark_df$mark+delay
    sample = paste(subjects,"-",visit,sep="") # putting subject and visit together
    return(sample)
}
# identify samples 2 visits after the dynamic mark
increase_dynamic_2m   = get_samples(increase_dynamic,delay=2)
stable_low_dynamic_2m = get_samples(stability_low_dynamic,delay=2)
stable_mod_dynamic_2m = get_samples(stability_moderate_dynamic,delay=2)

### CREATE DATA.FRAME FOR ANALYSIS
# setting variable type
info$Batch   = as.factor(info$Batch) 
info$Subject = as.factor(info$Subject) 

# selecting samples
# all dynamics
info_all = info[c(increase_dynamic_2m,
                  stable_low_dynamic_2m,
                  stable_mod_dynamic_2m),]
info_all$dynamic  = c(rep("increase",length(increase_dynamic_2m)),
                      rep("stable_low",length(stable_low_dynamic_2m)),
                      rep("stable_mod",length(stable_mod_dynamic_2m)))
info_all$dynamic = as.factor(info_all$dynamic)
# set1: stable low + increase
info_set1 = info[c(increase_dynamic_2m,
                   stable_low_dynamic_2m),]
info_set1$dynamic  = c(rep("increase",length(increase_dynamic_2m)),
                       rep("stable_low",length(stable_low_dynamic_2m)))
info_set1$dynamic = as.factor(info_set1$dynamic)
# set2: stable low + stable moderate
info_set2 = info[c(stable_mod_dynamic_2m,
                   stable_low_dynamic_2m),]
info_set2$dynamic  = c(rep("stable_mod",length(stable_mod_dynamic_2m)),
                       rep("stable_low",length(stable_low_dynamic_2m)))
info_set2$dynamic = as.factor(info_set2$dynamic)

# get voom processed data!
design_matrix = ~ ACE*PSS_2m_prior + (1|Batch) + (1|Subject)
unfiltered_design_matrix = ~ ACE*PSS_2m_prior + (1|Batch) + (1|Subject)

# identifying samples with PSS avaiable two months before!
complete_PSS_2m_samples= rownames(na.omit(info[,c("Sample","PSS_2m_prior")]))

voom_expr_unf = voomWithDreamWeights(counts = sncRNA_DGE[,complete_PSS_2m_samples], 
                                    formula = unfiltered_design_matrix, 
                                    data = info[complete_PSS_2m_samples,], 
                                    lib.size = info[complete_PSS_2m_samples,]$sequence_count)

voom_expr_all = voomWithDreamWeights(counts = sncRNA_DGE[,info_all$Sample], 
                                     formula = design_matrix, 
                                     data = info_all, 
                                     lib.size = info_all$sequence_count)

voom_expr_set1 = voomWithDreamWeights(counts = sncRNA_DGE[,info_set1$Sample], 
                                      formula = design_matrix, 
                                      data = info_set1, 
                                      lib.size = info_set1$sequence_count)

voom_expr_set2 = voomWithDreamWeights(counts = sncRNA_DGE[,info_set2$Sample], 
                                      formula = design_matrix, 
                                      data = info_set2, 
                                      lib.size = info_set2$sequence_count)

### WGCNA
# recommended setup in WGCNA documentation
options(stringsAsFactors = FALSE);
enableWGCNAThreads(nThreads = 12)

# Set of soft-thresholding power to check
powers = c(seq(4,10,by=1), seq(12,20, by=2));

# function to remove species notation in sncRNA name
removing_species_notation = function(data){
    names_sncRNA = gsub("hsa-","",colnames(data))
    names_sncRNA = gsub("Homo_sapiens_","",names_sncRNA)
    colnames(data) = names_sncRNA
    return(data)
}

# --------------------------------------------------------- #

## Set 1: Stable low + Increase
data_for_WGCNA = t(voom_expr_set1$E)
data_for_WGCNA = removing_species_notation(data_for_WGCNA)
cat("Samples:",nrow(data_for_WGCNA)) # 38

powerTables = list(data = pickSoftThreshold(data_for_WGCNA, powerVector=powers,verbose = F)[[2]]);
collectGarbage(); # select power 9

# WGCNA implementation
net = blockwiseModules(data_for_WGCNA, power = 9, networkType = "signed",
                        TOMType = "signed", minModuleSize = 19, # minModule = 10% of total
                        reassignThreshold = 0, mergeCutHeight = 0.25,
                        numericLabels = TRUE, pamRespectsDendro = FALSE,
                        verbose = 1)
# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
names(mergedColors) = names(net$colors)
# Plot the dendrogram and the module colors underneath
options(repr.plot.width = 6, repr.plot.height = 6)
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    main = "Set 1: Stable low + Increase \n(Power = 9; SFT_R2 = 0.868; mean.k = 2.78)",
                    "Cluster colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
# Calculate MEs with color labels
MEs0 = moduleEigengenes(data_for_WGCNA, mergedColors)$eigengenes
MEs = orderMEs(MEs0)
# cluster membership
datKME=signedKME(data_for_WGCNA, MEs, outputColumnName="MM_")
# Get cluster nodes
brown_1 = names(mergedColors[mergedColors == "brown"])
blue_1 = names(mergedColors[mergedColors == "blue"])
turquoise_1 = names(mergedColors[mergedColors == "turquoise"])
yellow_1 = names(mergedColors[mergedColors == "yellow"])

# create data.frame for statistical modeling
cluster_data =  data.frame(Subject = info_set1$Subject,
                           Batch = info_set1$Batch,
                           Visit = info_set1$Visit, 
                           PSS_2m_prior = info_set1$PSS_2m_prior,
                           dynamic  = info_set1$dynamic, 
                           ACE = info_set1$ACE,
                           Blue_activity = MEs$MEblue, 
                           Turquoise_activity = MEs$MEturquoise, 
                           Brown_activity = MEs$MEbrown, 
                           Yellow_activity = MEs$MEyellow)
rownames(cluster_data) = rownames(data_for_WGCNA)

# checking co-expression modules  
# ex. BLUE CLUSTER MODEL

# looking at the top 20 highest cluster membership nodes (kME) within each cluster
cat("BLUE\n")
datKME %>% arrange(desc(MM_blue)) %>% rownames() %>% head(,n=20) # depends on dplyr package

model_blue = lmer("Blue_activity ~ PSS_2m_prior + ACE + ACE:PSS_2m_prior + (1|Subject) + (1|Batch)",
                  cluster_data)

options(repr.plot.width = 10, repr.plot.height = 5)
# evaluating residual plot (DHARMa package)
temp = DHARMa::simulateResiduals(fittedModel = model_blue, plot = TRUE)
# checking effects
summary(model_blue)$coefficients # quick look at t-values for promissing significant effect 

# tab_model() function for mixed model result table (sjPlot package)
sjPlot::tab_model(model_blue, show.re.var = FALSE,digits = 5,p.threshold = c(0.05, 0.01, 0.001))

# --------------------------------------------------------- #

## Set 1: Stable low + Increase
data_for_WGCNA = t(voom_expr_set2$E)
data_for_WGCNA = removing_species_notation(data_for_WGCNA)
cat("Samples:",nrow(data_for_WGCNA)) # 53

powerTables = list(data = pickSoftThreshold(data_for_WGCNA, powerVector=powers,verbose = F)[[2]]);
collectGarbage(); # select power 8

# WGCNA implementation
net = blockwiseModules(data_for_WGCNA, power = 8, networkType = "signed",
                        TOMType = "signed", minModuleSize = 19, # minModule = 10% of total
                        reassignThreshold = 0, mergeCutHeight = 0.25,
                        numericLabels = TRUE, pamRespectsDendro = FALSE,
                        verbose = 1)
# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
names(mergedColors) = names(net$colors)
# Plot the dendrogram and the module colors underneath
options(repr.plot.width = 6, repr.plot.height = 6)
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    main = "Set 2: Stable low + Stable moderate \n(Power = 8; SFT_R2 = 0.756; mean.k = 3.72)",
                    "Cluster colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
# Calculate MEs with color labels
MEs0 = moduleEigengenes(data_for_WGCNA, mergedColors)$eigengenes
MEs = orderMEs(MEs0)
# cluster membership
datKME=signedKME(data_for_WGCNA, MEs, outputColumnName="MM_")
# Get cluster nodes
brown_2 = names(mergedColors[mergedColors == "brown"])
blue_2 = names(mergedColors[mergedColors == "blue"])
turquoise_2 = names(mergedColors[mergedColors == "turquoise"])
yellow_2 = names(mergedColors[mergedColors == "yellow"])

# create data.frame for statistical modeling
cluster_data =  data.frame(Subject = info_set2$Subject,
                           Batch = info_set2$Batch,
                           Visit = info_set2$Visit, 
                           PSS_2m_prior = info_set2$PSS_2m_prior,
                           dynamic  = info_set2$dynamic, 
                           ACE = info_set2$ACE,
                           Blue_activity = MEs$MEblue, 
                           Turquoise_activity = MEs$MEturquoise, 
                           Brown_activity = MEs$MEbrown, 
                           Yellow_activity = MEs$MEyellow)
rownames(cluster_data) = rownames(data_for_WGCNA)

# check co-expression modules

# --------------------------------------------------------- #
# Set 1 and Set 2 identified comparable co-expression communities...
# why not pool them together

## All: Stable low + Stable moderate + Increase

data_for_WGCNA = t(voom_expr_all$E)
data_for_WGCNA = removing_species_notation(data_for_WGCNA)
cat("Samples:",nrow(data_for_WGCNA)) # 63

powerTables = list(data = pickSoftThreshold(data_for_WGCNA, powerVector=powers,verbose = F)[[2]]);
collectGarbage(); # select power 8

# WGCNA implementation
net = blockwiseModules(data_for_WGCNA, power = 8, networkType = "signed",
                        TOMType = "signed", minModuleSize = 19, # minModule = 10% of total
                        reassignThreshold = 0, mergeCutHeight = 0.25,
                        numericLabels = TRUE, pamRespectsDendro = FALSE,
                        verbose = 1)
# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
names(mergedColors) = names(net$colors)
# Plot the dendrogram and the module colors underneath
options(repr.plot.width = 6, repr.plot.height = 6)
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    main = "Set 3: ALL \n(Power = 8; SFT_R2 = 0.757; mean.k = 3.68)",
                    "Cluster colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
# Calculate MEs with color labels
MEs0 = moduleEigengenes(data_for_WGCNA, mergedColors)$eigengenes
MEs = orderMEs(MEs0)
# cluster membership
datKME=signedKME(data_for_WGCNA, MEs, outputColumnName="MM_")
# Get cluster nodes
brown_3 = names(mergedColors[mergedColors == "brown"])
blue_3 = names(mergedColors[mergedColors == "blue"])
turquoise_3 = names(mergedColors[mergedColors == "turquoise"])
yellow_3 = names(mergedColors[mergedColors == "yellow"])

# create data.frame for statistical modeling
cluster_data =  data.frame(Subject = info_all$Subject,
                           Batch = info_all$Batch,
                           Visit = info_all$Visit, 
                           PSS_2m_prior = info_all$PSS_2m_prior,
                           dynamic  = info_all$dynamic, 
                           ACE = info_all$ACE,
                           Blue_activity = MEs$MEblue, 
                           Turquoise_activity = MEs$MEturquoise, 
                           Brown_activity = MEs$MEbrown, 
                           Yellow_activity = MEs$MEyellow)
rownames(cluster_data) = rownames(data_for_WGCNA)

# --------------------------------------------------------- #

## Unfiltered WGCNA 
data_for_WGCNA = t(voom_expr_unf$E)
data_for_WGCNA = removing_species_notation(data_for_WGCNA)
cat("Samples:",nrow(data_for_WGCNA)) # 175, so highest sample size and statistical power

powerTables = list(data = pickSoftThreshold(data_for_WGCNA, powerVector=powers,verbose = F)[[2]]);
collectGarbage(); # select power 6

# WGCNA implementation
net = blockwiseModules(data_for_WGCNA, power = 6, networkType = "signed",
                        TOMType = "signed", minModuleSize = 19, # minModule = 10% of total
                        reassignThreshold = 0, mergeCutHeight = 0.25,
                        numericLabels = TRUE, pamRespectsDendro = FALSE,
                        verbose = 1)
# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
to_brown  = which(mergedColors == "yellow")
to_yellow = which(mergedColors == "brown")
mergedColors[to_brown] = "brown"
mergedColors[to_yellow] = "yellow"

names(mergedColors) = names(net$colors)
# Plot the dendrogram and the module colors underneath
options(repr.plot.width = 6, repr.plot.height = 6)
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    main = "All samples collected \n(Power = 6; SFT_R2 = 0.815; median.k = 4.73)",
                    "Cluster colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
# Calculate MEs with color labels
MEs0 = moduleEigengenes(data_for_WGCNA, mergedColors)$eigengenes
MEs = orderMEs(MEs0)
# cluster membership
datKME=signedKME(data_for_WGCNA, MEs, outputColumnName="MM_")
# Get cluster nodes
brown_4 = names(mergedColors[mergedColors == "brown"])
blue_4 = names(mergedColors[mergedColors == "blue"])
turquoise_4 = names(mergedColors[mergedColors == "turquoise"])
yellow_4 = names(mergedColors[mergedColors == "yellow"])

options(repr.plot.width = 5, repr.plot.height = 5)
diss1=1-TOMsimilarityFromExpr(data_for_WGCNA, power =6)
diag(diss1) = NA


# TOM PLOT
#pdf(file = "network_TOM_all_samples.pdf",width = 5,height = 5)
a = TOMplot(diss1^6, net$dendrograms[[1]], mergedColors,
        mergedColors,
        main = "Network TOM and co-expression clusters")
#dev.off()

# Function to get sncRNA type composition
print_prop = function(set){
    cat("Total size:",length(set),"\n")
    cat("miRNA: ",paste(length(grep("miR|let",set)),
               " (",length(grep("miR|let",set))/length(set),sep=""),")\n",sep="")
    cat("piRNA: ",paste(length(grep("piR",set)),
               " (",length(grep("piR",set))/length(set),sep=""),")\n",sep="")
    cat("tRNA: ",paste(length(grep("tR",set)),
               " (",length(grep("tR",set))/length(set),sep=""),")\n",sep="")
}
# example
print_prop(brown_4)
#Total size: 32 
#miRNA: 18 (0.5625)
#piRNA: 14 (0.4375)
#tRNA: 0 (0)

# PIE CHART WITH sncRNA TYPE COMPOSITION IN EACH CO-EXPERSSION MODULE
colrs =c('hotpink3','steelblue','gold') # miRNA; piRNA; tRNA
#pdf(file = "brow_cluster_composition.pdf",width = 5,height = 5)
pie(c(18,14),labels =  c("18 miRNAs","14 piRNAs"),
    radius = 0.5, sub = "Brown cluster",
    border = "white", col=c('hotpink3','steelblue'))
#dev.off()
#pdf(file = "blue_cluster_composition.pdf",width = 5,height = 5)
pie(c(3,37),labels =  c("3 miRNAs","37 piRNAs"),
    radius = 0.5, sub = "Blue cluster",
    border = "white", col=c('hotpink3','steelblue'))
#dev.off()
#pdf(file = "turquoise_cluster_composition.pdf",width = 5,height = 5)
pie(c(14,55,2),labels =  c("14 miRNAs","55 piRNAs","2 tRNAs"),
    radius = 0.5, sub = "Turquoise cluster",
    border = "white", col=c('hotpink3','steelblue','gold'))
#dev.off()
#pdf(file = "yellow_cluster_composition.pdf",width = 5,height = 5)
pie(c(5,23,8),labels =  c("5 miRNAs","23 piRNAs","8 tRNAs"),
    radius = 0.5, sub = "Yellow cluster",
    border = "white", col=c('hotpink3','steelblue','gold'))
#dev.off()


## PLOT CLUSTER REPRESENTATIVE sncRNAs
datKME$class = 0 
datKME[grep("let|miR",rownames(datKME)),"class"] = "miRNA"
datKME[grep("tR",rownames(datKME)),"class"]      = "tRNA"
datKME[grep("piR",rownames(datKME)),"class"]     = "piRNA"
datKME$sncRNA = rownames(datKME)

plot1 = ggplot(datKME[!(datKME$MM_turquoise < 0.8),], aes(y = forcats::fct_reorder(sncRNA,MM_turquoise), 
                                     x = MM_turquoise, fill = class)) + geom_col(alpha=0.9) + 
ylab("") + xlab("Membership") +
ggtitle("Turquoise cluster representative sncRNAs", "Ranking cluster membership") + theme_bw() + 
theme(plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 12),
        axis.title.x = element_text(size = 12),
        axis.text.y = element_text(size = 10.5),
        axis.text.x = element_text(size = 12),
        legend.position="top") +
        geom_vline(xintercept = 0.8,linetype="dashed") +
scale_x_continuous(limits = c(0, 1), breaks = seq(0,1,0.2)) +
scale_fill_manual(values=colrs)

#pdf(file = "ranking_representative_sncRNAs_turquoise.pdf",width = 6,height = 7)
plot1
#dev.off()

# create data.frame for statistical modeling
cluster_data =  data.frame(Subject = info[complete_PSS_2m_samples,]$Subject,
                           Batch = info[complete_PSS_2m_samples,]$Batch,
                           Visit = info[complete_PSS_2m_samples,]$Visit,
                           PSS = info[complete_PSS_2m_samples,]$PSS_Score,
                           PSS_delta_acute = info[complete_PSS_2m_samples,]$PSS_delta_acute,
                           PSS_2m_prior = info[complete_PSS_2m_samples,]$PSS_2m_prior,
                           PSS_delta_2m_delayed = info[complete_PSS_2m_samples,]$PSS_delta_2m_delayed,
                           ACE = info[complete_PSS_2m_samples,]$ACE,
                           Blue_activity = MEs$MEblue, 
                           Turquoise_activity = MEs$MEturquoise, 
                           Brown_activity = MEs$MEbrown, 
                           Yellow_activity = MEs$MEyellow)
rownames(cluster_data) = rownames(data_for_WGCNA)


## Turquoise cluster model

model_turquoise = lmer("Turquoise_activity ~ PSS_2m_prior*ACE  + (1|Subject) + (1|Batch)",
                  cluster_data)
options(repr.plot.width = 10, repr.plot.height = 5)
# evaluating residual plot (DHARMa package)
#pdf(file = "DHARMa_residual_simulation_linear_model.pdf",width = 10,height = 5)
temp = DHARMa::simulateResiduals(fittedModel = model_turquoise, plot = TRUE)
# dev.off()
# checking effects
summary(model_turquoise)$coefficients

sjPlot::tab_model(model_turquoise, show.re.var = FALSE, show.stat = TRUE,show.aic = TRUE,
                  digits = 4,p.threshold = c(0.05, 0.01, 0.001))

# addressing heteroskedasticity
model_turquoise_corrected = lmer("Turquoise_activity ~ I(PSS_2m_prior^2) + I(ACE^2) + PSS_2m_prior:ACE + (1|Subject) + (1|Batch)",
                                  cluster_data[])
options(repr.plot.width = 10, repr.plot.height = 5)
# evaluating residual plot (DHARMa package)
# pdf(file = "DHARMa_residual_simulation_quadratic_model.pdf",width = 10,height = 5)
temp = DHARMa::simulateResiduals(fittedModel = model_turquoise_corrected, plot = TRUE)
# dev.off()
# checking effects
summary(model_turquoise_corrected)$coefficients

sjPlot::tab_model(model_turquoise_corrected, show.stat = TRUE,show.aic = TRUE,
show.re.var = FALSE,digits = 4,p.threshold = c(0.05, 0.01, 0.001))

performance::compare_performance(model_turquoise,model_turquoise_corrected) # comparing models using performance package

options(repr.plot.width = 4, repr.plot.height = 4,repr.plot.res = 100)

## PLOTING EFFECT
# pdf(file = "turquoise_ACE_effect.pdf",width = 4,height = 4)
plot(effects::Effect("ACE", model_turquoise_corrected,
                     xlevels=list(ACE=seq(0, 9, 1))), 
     main = expression("ACE effect plot"), lwd=2,
     ylab = "Turquoise cluster activity",
     xlab = expression("ACE"),grid=TRUE)
# dev.off()

# pdf(file = "turquoise_PSS_effect.pdf",width = 4,height = 4)
plot(effects::Effect("PSS_2m_prior", model_turquoise_corrected,
                     xlevels=list(PSS_2m_prior=seq(0, 30, 1))), 
     main = expression("PSS"[t-2]~"effect plot"), lwd=2,
     ylab = "Turquoise cluster activity",
     xlab = expression("PSS"[t-2]),grid=TRUE)
# dev.off()


