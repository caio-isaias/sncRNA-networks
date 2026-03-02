### LOAD PACKAGES
library(WGCNA)
library(edgeR)
library(variancePartition)

### LOAD DATA
path = "analysis_manuscript_DREAM_WGCNA/self_contained_analysis/data"

miRNA_DGE = readRDS(paste(path,"/DGEList_miRNA_filtered_norm.RDS",sep=""))
piRNA_DGE = readRDS(paste(path,"/DGEList_piRNA_filtered_norm.RDS",sep=""))
tRNA_DGE  = readRDS(paste(path,"/DGEList_tRNA_filtered_norm.RDS",sep=""))

dynamic_candidates  = readRDS(paste(path,"/dynamic_cadidate_all_subjects_DC_DE.RDS",sep=""))

info = read.csv(paste(path,"/data_for_modeling.csv",sep=""),
                row.names = read.csv(paste(path,"/data_for_modeling.csv",sep=""))[,1])[,-1]

removing_species_notation = function(data){
    names_sncRNA = gsub("hsa-","",colnames(data))
    names_sncRNA = gsub("Homo_sapiens_","",names_sncRNA)
    colnames(data) = names_sncRNA
    return(data)
}

# filter out non-dynamic sncRNA
dynamic_miRNA = dynamic_candidates$miRNA_DE[(dynamic_candidates$miRNA_DE$dynamic_count_DE >= 7),]$miRNA
dynamic_piRNA = dynamic_candidates$piRNA_DE[(dynamic_candidates$piRNA_DE$dynamic_count_DE >= 7),]$piRNA
dynamic_tRNA = dynamic_candidates$tRNA_DE[(dynamic_candidates$tRNA_DE$dynamic_count_DE >= 7),]$tRNA

sncRNA_DGE = rbind(miRNA_DGE[dynamic_miRNA,rownames(info)],
                   piRNA_DGE[dynamic_piRNA,rownames(info)],
                   tRNA_DGE[dynamic_tRNA,rownames(info)])

# get the dynamic PSS reference visits
increase_dynamic = read.csv(paste(path,"/mark_increase_dynamic.csv",sep=""))[,-1]
stability_low_dynamic = read.csv(paste(path,"/mark_stability_low_dynamic.csv",sep=""))[,-1]
stability_moderate_dynamic = read.csv(paste(path,"/mark_stability_moderate_dynamic.csv",sep=""))[,-1]

# function to get samples from PSS dynamic classification data.frame  
get_samples <- function(mark_df, delay = 0){
    subjects = gsub("50","",mark_df$Record_ID) # remove "50", from the Record_ID
    visit = mark_df$mark+delay
    sample = paste(subjects,"-",visit,sep="") # putting subject and visit together
    return(sample)
}
# identify samples 2 visits after the dynamic
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

# set1: stable low 
info_low = info[c(stable_low_dynamic_2m),]
# set2:  stable moderate
info_mod = info[c(stable_mod_dynamic_2m),]
# set3: increase
info_incres = info[c(increase_dynamic_2m),]

# get voom processed data!
design_matrix_null = ~  (1|Batch) + (1|Subject)
design_matrix = ~ ACE*PSS_2m_prior + (1|Batch) + (1|Subject)
design_matrix2 = ~ ACE*PSS_2m_prior + (1|Batch) 

# identifying samples with PSS avaiable two months before!
complete_PSS_2m_samples= rownames(na.omit(info[,c("Sample","PSS_2m_prior")]))

voom_expr_unf = voomWithDreamWeights(counts = sncRNA_DGE[,complete_PSS_2m_samples], 
                                    formula = design_matrix, 
                                    data = info[complete_PSS_2m_samples,], 
                                    lib.size = info[complete_PSS_2m_samples,]$sequence_count)

voom_expr_all = voomWithDreamWeights(counts = sncRNA_DGE[,info_all$Sample], 
                                     formula = design_matrix, 
                                     data = info_all, 
                                     lib.size = info_all$sequence_count)

voom_expr_set1 = voomWithDreamWeights(counts = sncRNA_DGE[,info_low$Sample], 
                                      formula = design_matrix, 
                                      data = info_low, 
                                      lib.size = info_low$sequence_count)

voom_expr_set2 = voomWithDreamWeights(counts = sncRNA_DGE[,info_mod$Sample], 
                                      formula = design_matrix, 
                                      data = info_mod, 
                                      lib.size = info_mod$sequence_count)
voom_expr_set3 = voomWithDreamWeights(counts = sncRNA_DGE[,info_incres$Sample], 
                                      formula = design_matrix2, 
                                      data = info_incres, 
                                      lib.size = info_incres$sequence_count)

# identify samples with high or low ACE score
ACE_low = rownames(info[info$ACE<=1,])
ACE_high = rownames(info[info$ACE>1,])

voom_expr_ACE_low = voomWithDreamWeights(counts = sncRNA_DGE[,ACE_low], 
                                     formula = design_matrix_null, 
                                     data = info[ACE_low,], 
                                     lib.size = info[ACE_low,]$sequence_count)

voom_expr_ACE_high = voomWithDreamWeights(counts = sncRNA_DGE[,ACE_high], 
                                     formula = design_matrix_null, 
                                     data = info[ACE_high,], 
                                     lib.size = info[ACE_high,]$sequence_count)

### REFERENCE WGCNA (LOW STABLE STRESS)
wgcna_low = t(voom_expr_set1$E) # stable low
wgcna_low = removing_species_notation(wgcna_low)
cat("Sample:",nrow(wgcna_low))

# Set of soft-thresholding power to check
powers = c(seq(4,10,by=1), seq(12,20, by=2));
powerTables = list(data = pickSoftThreshold(wgcna_low, powerVector=powers,verbose = F)[[2]]);
collectGarbage(); # select power 9

## WGCNA
net = blockwiseModules(wgcna_low, power = 14, networkType = "signed",
                        TOMType = "signed", minModuleSize = 19, # minModule = 10% of total
                        reassignThreshold = 0, mergeCutHeight = 0.25,
                        numericLabels = TRUE, pamRespectsDendro = FALSE,
                        verbose = 1)
# Convert labels to colors for plotting
reference_colors = labels2colors(net$colors)
names(reference_colors) = names(net$colors)
multiColor  = list(Set1 = reference_colors)
# Plot the dendrogram and the module colors underneath
options(repr.plot.width = 6, repr.plot.height = 6)
plotDendroAndColors(net$dendrograms[[1]], reference_colors[net$blockGenes[[1]]],
                    main = "Stable dynamic with low stress \n(Power = 14; SFT_R2 = 0.840; median.k = 1.03)",
                    "Cluster colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
# Calculate MEs with color labels
MEs0 = moduleEigengenes(wgcna_low, reference_colors)$eigengenes
MEs = orderMEs(MEs0)
# cluster membership
datKME=signedKME(wgcna_low, MEs, outputColumnName="MM_")
# Get cluster nodes
brown = names(reference_colors[reference_colors == "brown"])
blue = names(reference_colors[reference_colors == "blue"])
turquoise = names(reference_colors[reference_colors == "turquoise"])
yellow = names(reference_colors[reference_colors == "yellow"])

## TOM PLOT
diss1=1-TOMsimilarityFromExpr(wgcna_low, power =14)
diag(diss1) = NA
#pdf(file = "reference_network_TOM.pdf",width = 5,height = 5)
a = TOMplot(diss1^14, net$dendrograms[[1]], reference_colors,
            main = "Network TOM and co-expression clusters")
#dev.off()

### PRESERVATION ANALYSIS
wgcna_all = t(voom_expr_all$E)
wgcna_all = removing_species_notation(wgcna_all)
wgcna_mod = t(voom_expr_set2$E) # stable moderate
wgcna_mod = removing_species_notation(wgcna_mod)
wgcna_inc = t(voom_expr_set3$E) # increase
wgcna_inc = removing_species_notation(wgcna_inc)

stable_multidata   = multiData(Set1 = wgcna_low, Set2 = wgcna_mod)
increase_multidata = multiData(Set1 = wgcna_low, Set2 = wgcna_inc)

# module preservation
preservation_in_stable =   modulePreservation(stable_multidata, 
                                              multiColor = multiColor, 
                                              dataIsExpr = TRUE,
                                              nPermutations = 100,
                                              networkType = "signed", 
                                              referenceNetworks = 1)

preservation_in_increase =   modulePreservation(increase_multidata, 
                                                multiColor = multiColor, 
                                                dataIsExpr = TRUE,
                                                nPermutations = 100,
                                                networkType = "signed", 
                                                referenceNetworks = 1)

# put preservation result in a data.set
result_stable = preservation_in_stable$preservation$Z$ref.Set1$inColumnsAlsoPresentIn.Set2[,c(2:4)]
result_stable = result_stable[-which(rownames(result_stable) %in% c("gold","grey")),]
colnames(result_stable) = c("Z-Summary","Z-Density","Z-Connectivity")
result_stable$module = rownames(result_stable)
result_stable_melted = reshape2::melt(result_stable, id.vars = "module")

result_increase = preservation_in_increase$preservation$Z$ref.Set1$inColumnsAlsoPresentIn.Set2[,c(2:4)]
result_increase = result_increase[-which(rownames(result_increase) %in% c("gold","grey")),]
colnames(result_increase) = c("Z-Summary","Z-Density","Z-Connectivity")
result_increase$module = rownames(result_increase)
result_increase_melted = reshape2::melt(result_increase, id.vars = "module")

# combine
result_melted = rbind(result_stable_melted,result_increase_melted)
result_melted$comparison = c(rep("Compared to stable moderate-stress",nrow(result_stable_melted)),
                             rep("Compared to increased-stress",nrow(result_increase_melted)))

### PLOT PRESERVATION RESULTS 
library(ggplot2)
# pdf(file = "module_preservation_test.pdf",width = 6.5,height = 6)
ggplot(result_melted, aes(x = variable, y = value, color = module)) +
  geom_point(size = 3) +
  geom_hline(yintercept = 2, color = "red", linetype = "dashed") +
  geom_hline(yintercept = 10, color = "darkgreen", linetype = "dashed") +
  theme_bw() + xlab("") +
  theme(axis.text.x = element_text(color="black",face = "bold"),
        axis.text.y = element_text(face = "bold"),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_color_manual(name = 'Module', values =c("blue",'brown','turquoise',"gold"), 
                            labels = c('blue','brown','turquoise',"yellow")) +
  labs(title = "Module Preservation test",
       subtitle = "Stable low-stress as reference network",
       y = "Z-Statistic") +
    facet_grid(~comparison)
# dev.off()
