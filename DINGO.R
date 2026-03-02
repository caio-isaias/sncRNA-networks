### LOAD PACKAGES
library("iDINGO")
library("WGCNA")
library("edgeR"")
library("igraph")
library("ggplot2")

### LOAD DATA
# path to files
path_nickole = "/media/caio/T7_Shield/PhD/BEPE/Backup/Tracy_Server_Backup/Nickole_Study_2024/data_preprocessing/"
miRNA_DGE = readRDS(paste(path_nickole,"DGEList_miRNA_filtered_norm.RDS",sep="")); 
piRNA_DGE = readRDS(paste(path_nickole,"DGEList_piRNA_filtered_norm.RDS",sep="")); 
tRNA_DGE  = readRDS(paste(path_nickole,"DGEList_tRNA_filtered_norm.RDS",sep=""));
# path to sncRNA dynamic classification results
dynamic_candidates  = readRDS(paste(path_nickole,"dynamic_cadidate_all_subjects_DC_DE.RDS",sep=""));

info_to_modeling = read.csv("data_for_modeling.csv", 
                            row.names = read.csv("data_for_modeling.csv")[,1])[,-1]

# filter out non-dynamic sncRNA
dynamic_miRNA = dynamic_candidates$miRNA_DE[(dynamic_candidates$miRNA_DE$dynamic_count_DE >= 7),]$miRNA
dynamic_piRNA = dynamic_candidates$piRNA_DE[(dynamic_candidates$piRNA_DE$dynamic_count_DE >= 7),]$piRNA
dynamic_tRNA = dynamic_candidates$tRNA_DE[(dynamic_candidates$tRNA_DE$dynamic_count_DE >= 7),]$tRNA

sncRNA_DGE = rbind(miRNA_DGE[dynamic_miRNA,rownames(info_to_modeling)],
                   piRNA_DGE[dynamic_piRNA,rownames(info_to_modeling)],
                   tRNA_DGE[dynamic_tRNA,rownames(info_to_modeling)])

path = "/home/caio/Desktop/workfolder/BEPE/analysis_manuscript_DREAM_WGCNA/self_contained_analysis/data"
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

# reference samples for ACE's.
ACE_samples = c("006-4","007-8","014-4","016-6","018-7","019-5","021-9","023-9",
                "025-9","037-5","042-4","043-6","046-4","055-4","068-4","072-4","074-7",
                "075-9","080-4","084-5","085-4","086-7","088-9","104-9","107-5","112-4",
                "121-4")

# remove subjects with high ACE (ACE = 9  and ACE = 7)
stable_mod_dynamic_2m = stable_mod_dynamic_2m[-grep("007|107",stable_mod_dynamic_2m)]

info_PSS = info_to_modeling[c(stable_low_dynamic_2m,
                              stable_mod_dynamic_2m),]
PSS_dynamic = as.factor(c(rep(0,length(stable_low_dynamic_2m)),
                          rep(1,length(stable_mod_dynamic_2m)))) # changes if we're comparing low vs. increased
ACE_factor_group = as.factor(as.numeric(info_to_modeling[ACE_samples,"ACE"]>1))

## initialize voom
# ACE
design_matrix_ACE = ~ ACE + (1|Batch) + (1|Subject)
voom_expr = variancePartition::voomWithDreamWeights(sncRNA_DGE, 
                                                    formula = design_matrix_ACE, 
                                                    info_to_modeling, 
                                                    lib.size = info_to_modeling$sequence_count)
data_for_dingo = t(voom_expr$E)
# PSS
design_matrix_PSS = ~ (1|Batch) + (1|Subject)
voom_expr = variancePartition::voomWithDreamWeights(sncRNA_DGE, 
                                                    formula = design_matrix_PSS, 
                                                    info_to_modeling, 
                                                    lib.size = info_to_modeling$sequence_count)
data_for_dingo2 = t(voom_expr$E)

### DINGO - PSS (LOW vs. MODERATE)
PSS_samples=c(stable_low_dynamic_2m,stable_mod_dynamic_2m)
# takes a lot of time to run this... (22h for 16gb ram)
dingo_result_voom2 = dingo(dat = data_for_dingo2[PSS_samples,],PSS_dynamic,  
                           cores= 12, B = 100,verbose=T)
saveRDS(dingo_result_voom2,"dingo_result_all_samples_PSS_stable.RDS")
#dingo_result_voom2 = readRDS("dingo_result_all_samples_PSS_stable.RDS")

# renaming sncRNAs
dingo_result_voom2$genepair$gene1 = gsub("hsa-","",dingo_result_voom2$genepair$gene1)
dingo_result_voom2$genepair$gene1 = gsub("Homo_sapiens_","",dingo_result_voom2$genepair$gene1)
dingo_result_voom2$genepair$gene2 = gsub("hsa-","",dingo_result_voom2$genepair$gene2)
dingo_result_voom2$genepair$gene2 = gsub("Homo_sapiens_","",dingo_result_voom2$genepair$gene2)


### DINGO - PSS (LOW vs. INCREASED)

PSS_samples=c(stable_low_dynamic_2m,increase_dynamic_2m)
PSS_dynamic = as.factor(c(rep(0,length(stable_low_dynamic_2m)),
                          rep(1,length(increase_dynamic_2m))))
dingo_result_voom2 = dingo(dat = data_for_dingo2[PSS_samples,],PSS_dynamic, cores= 12, B = 100,verbose=T)
saveRDS(dingo_result_voom2,"dingo_result_all_samples_PSS_increased.RDS")
dingo_result_increase = readRDS("dingo_result_all_samples_PSS_increase.RDS")

# prepare network
signf_id = dingo_result_increase$p.val<0.001
edges = data.frame(gene1 = dingo_result_increase$genepair$gene1[signf_id],
                   gene2 = dingo_result_increase$genepair$gene2[signf_id],
                   weight = dingo_result_increase$diff.score[signf_id])
diff_network = graph_from_data_frame(edges, directed = F) # create graph with iGraph function
# rename sncRNAs
miRNA_names = gsub("hsa-","",miRNA_DGE$genes$genes)
piRNA_names = gsub("hsa-","",piRNA_DGE$genes$genes)
tRNA_names = gsub("Homo_sapiens_","",tRNA_DGE$genes$genes)
# rename sncRNA in the network
V(diff_network)$name = gsub("hsa-","",V(diff_network)$name)
V(diff_network)$name = gsub("Homo_sapiens_","",V(diff_network)$name)
class = rep(0,length(V(diff_network)$name))
class[V(diff_network)$name %in% miRNA_names] = "miRNA"
class[V(diff_network)$name %in% piRNA_names] = "piRNA"
class[V(diff_network)$name %in% tRNA_names] = "tRNA"

diff_network <- set_vertex_attr(diff_network, "class", value = as.factor(class))
colrs <- c("hotpink3","steelblue", "gold")
V(diff_network)$color <- colrs[V(diff_network)$class]
# Compute node degrees and use that to set node size:
deg <- degree(diff_network, mode="all")
V(diff_network)$size <- ((max(deg)+deg)/max(deg))*4
# Set edge width based on weight:
E(diff_network)$width <- abs(E(diff_network)$weight)*0.9

options(repr.plot.width = 10, repr.plot.height = 10)

# pdf(file = "differential_network_DINGO_PSS_stable.pdf",width = 7,height = 7)
plot(diff_network, layout=layout.circle, 
     edge.color=ifelse(E(diff_network)$weight > 0, "blue","firebrick2"),
     asp = 0,vertex.label=NA,cex.main=2)
title("Differential Network \n(Stable low-stress vs. Increased-stress; p<0.001)",cex.main=1.5)
legend(x="bottomleft", c("miRNA","piRNA","tRNA"), pch=21,
       col="#777777", pt.bg=colrs, pt.cex=2, cex=1, bty="n", ncol=1)
# dev.off()

total_degree = degree(diff_network, mode="all")
degree_df = data.frame(sncRNA = V(diff_network)$name,
                       degree_total = total_degree,
                       negative = 0,
                       positive = 0,
                       class = class)

# positive & negative edges
positive_edges = edges[edges$weight > 0,]
negative_edges = edges[edges$weight < 0,]

for(gene in V(diff_network)$name){
    positive = sum(positive_edges[which(positive_edges$gene1 == gene | positive_edges$gene2 == gene),]$weight)
    negative = sum(negative_edges[which(negative_edges$gene1 == gene | negative_edges$gene2 == gene),]$weight)
    
    degree_df[degree_df$sncRNA == gene,"positive"] = positive
    degree_df[degree_df$sncRNA == gene,"negative"] = negative
}

rename_edge = function(x){
    x$gene1 = gsub("hsa-","",x$gene1) 
    x$gene2 = gsub("hsa-","",x$gene2)
    x$gene1 = gsub("Homo_sapiens_","",x$gene1) 
    x$gene2 = gsub("Homo_sapiens_","",x$gene2)
    return(x)
}
positive_edges = rename_edge(positive_edges)
negative_edges = rename_edge(negative_edges)

options(repr.plot.width = 6, repr.plot.height = 10)


plot1 = ggplot(degree_df[!(degree_df$degree_total < 5),], aes(y = forcats::fct_reorder(sncRNA,degree_total), 
                                     x = degree_total, fill = class)) + geom_col(alpha=0.9) + 
ylab("") + xlab("Degree") +
ggtitle("Differential network hub sncRNA", "Ranking most connected sncRNA") + theme_bw() + 
theme(plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 12),
        axis.title.x = element_text(size = 12),
        axis.text.y = element_text(size = 10.5),
        axis.text.x = element_text(size = 12),
        legend.position="top") +
        #geom_vline(xintercept = 0.9,linetype="dashed")
scale_x_continuous(limits = c(0, 65), breaks = c(0,10,20,30,40,50,60,70)) +
scale_fill_manual(values=colrs)

plot2 = ggplot(degree_df[!(degree_df$positive < 10),], aes(y = forcats::fct_reorder(sncRNA,positive), 
                                     x = positive, fill = class)) + geom_col(alpha=0.9) + 
ylab("") + xlab("Weighted degree") +
ggtitle("sncRNA sum of positive edges", "Ranking sncRNA with the highest positive differential score") + theme_bw() + 
theme(plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 12),
        axis.title.x = element_text(size = 12),
        axis.text.y = element_text(size = 10.5),
        axis.text.x = element_text(size = 12),
        legend.position="none") +
        #geom_vline(xintercept = 0.9,linetype="dashed")
scale_x_continuous(limits = c(0, 60), breaks = c(0,10,20,30,40,50,60)) +
scale_fill_manual(values=colrs)

plot3 = ggplot(degree_df[!(degree_df$negative > -10),], aes(y = forcats::fct_reorder(sncRNA,negative, .desc = TRUE), 
                                     x = negative, fill = class)) + geom_col(alpha=0.9) + 
ylab("") + xlab("Weighted degree") +
ggtitle("sncRNA sum of negative edges", "Ranking sncRNA with the highest negative differential score") + theme_bw() + 
theme(plot.title = element_text(hjust = 0.5, size = 15, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 12),
        axis.title.x = element_text(size = 12),
        axis.text.y = element_text(size = 10.5),
        axis.text.x = element_text(size = 12),
        legend.position="none") +
        #geom_vline(xintercept = 0.9,linetype="dashed")
scale_x_continuous(limits = c(0,-70), breaks = c(0,-10,-20,-30,-40,-50,-60,-70)) +
scale_x_reverse() +
scale_fill_manual(values=colrs)

gridExtra::grid.arrange(plot1, gridExtra::arrangeGrob(plot2, plot3, heights=c(1,1)), ncol = 2)

### FUNCTION TO PLOT SUBGRAPH FOR A SINGLE sncRNA
plot_subgraph <- function(sncRNA,positive = TRUE, 
                          node_scale = 4, edge_scale = 0.9,
                          dist = 1.8){
    if(positive){
        table = positive_edges[(positive_edges$gene1 == sncRNA)|(positive_edges$gene2 == sncRNA),]
        edge_color = "blue"
        title= paste(sncRNA," - positive score",sep="")
    }else{
        table = negative_edges[(negative_edges$gene1 == sncRNA)|(negative_edges$gene2 == sncRNA),]
        edge_color = "red"
        title= paste(sncRNA," - negative score",sep="")
    }
    diff_subnetwork = graph_from_data_frame(table, directed = F)
    
    n = vcount(diff_subnetwork)
    class = rep(0,n)
    class[V(diff_subnetwork)$name %in% miRNA_names] = "miRNA"
    class[V(diff_subnetwork)$name %in% piRNA_names] = "piRNA"
    class[V(diff_subnetwork)$name %in% tRNA_names] = "tRNA"

    diff_subnetwork <- set_vertex_attr(diff_subnetwork, "class", value = as.factor(class))
    colrs <- c("hotpink3","steelblue", "gold")
    V(diff_subnetwork)$color <- colrs[V(diff_subnetwork)$class]

    # Compute node degrees and use that to set node size:
    deg <- degree(diff_subnetwork, mode="all")
    V(diff_subnetwork)$size <- ((max(deg)+deg)/max(deg))*node_scale

    # Set edge width based on weight:
    E(diff_subnetwork)$width <- abs(E(diff_subnetwork)$weight)*edge_scale 
    radian.rescale <- function(x, start=0, direction=1){
      c.rotate <- function(x) (x + start) %% (2 * pi) * direction
      c.rotate(scales::rescale(x, c(0, 2 * pi), range(x)))
    }
    lab.locs <- radian.rescale(x=1:n, direction=-1, start=0)
    plot(diff_subnetwork, layout=layout.circle, 
     edge.color=edge_color,
     asp = 0,cex.main=2, 
     vertex.label.dist = dist,
     vertex.label.degree = lab.locs,
     vertex.label.font = 2,
     vertex.label.color = "black")
    title(title,cex.main=1.5)
    legend = FALSE
}

## EXAMPLE:
# pdf(file = "DINGO_subgraph_miR_449a_negative_score.pdf",width = 6,height = 6)
plot_subgraph("miR-449a",positive = FALSE, dist=2, node_scale =8)
# dev.off()
# pdf(file = "DINGO_subgraph_miR_449a_positive_score.pdf",width = 6,height = 6)
plot_subgraph("miR-449a",positive = TRUE, dist=, node_scale =8)
# dev.off()

