## Figure 4
## code to plot mirrorplots for nia-aa-a-score in dentate gyrus
## this will use promoter associated protein coding ENS-IDs only
library(tidyverse)
library(hudson)
library(RnBeads)
library(grid)
library(ggpubr)

setwd("~/90plus")
baseDir <- getwd()
savedir <- paste0("~/90plus/")

## get annotation from rnbeads
anno_promoters <- rnb.annotation2data.frame(rnb.get.annotation("promoters"))

## specify region 
all_regions <- c("DG")

## functions

## function to make a manhattan mirror plot, code adapted from R package gmirror.
## https://rdrr.io/github/anastasia-lucas/hudson/src/R/gmirror.R 
adapted_emirror <- function(top, bottom, tline, bline, log_10=TRUE, yaxis, opacity=0.5, chroms = c(1:22),
                            toptitle=NULL, bottomtitle=NULL, annotate_var, annotate_p,
                            highlight_var, highlight_p, highlighter="skyblue", color1="#AAAAAA", 
                            color2="#4D4D4D", groupcolors, rotatelabels=FALSE, labelangle = 90, 
                            freey=FALSE, background="white", grpblocks=FALSE, 
                            file="emirror", type="png", hgtratio=2, hgt=7, wi=12, res=300){
  
  ###### Combine data ######
  topn <- names(top)
  bottomn <- names(bottom)
  top$Location <- "Top"
  bottom$Location <- "Bottom"
  
  ### File format check ###
  if(!identical(topn, bottomn)){stop("Please ensure both inputs have the same metadata columns.")}
  d <- as.data.frame(rbind(top, bottom))
  
  ####### Info for y-axis ######
  if(log_10==TRUE){
    d$pval <- -log10(d$pvalue)
    yaxislab1 <- expression(paste("-log"[10], "(p)", sep=""))
    yaxislab2 <- expression(paste("-log"[10], "(p)", sep=""))
    if(!missing(tline)) {tredline <- -log10(tline)}
    if(!missing(bline)) {bredline <- -log10(bline)}
  } else {
    d$pval <- d$pvalue
    yaxislab1 <- yaxis[1]
    yaxislab2 <- yaxis[2]
    if(!missing(tline)) {tredline <- tline}
    if(!missing(bline)) {bredline <- bline}
  }
  yaxismax1 <- ifelse(freey==FALSE, max(d$pval[which(d$pval< Inf)]), max(d$pval[which(d$pval< Inf) & d$Location=="Top"]))
  yaxismax2 <- ifelse(freey==FALSE, max(d$pval[which(d$pval< Inf)]), max(d$pval[which(d$pval< Inf) & d$Location=="Bottom"]))
  yaxismin1 <- ifelse(freey==FALSE, 0, min(d$pval[d$Location=="Top"]))
  yaxismin2 <- ifelse(freey==FALSE, 0, min(d$pval[d$Location=="Bottom"]))
  
  ###### Save to merge later ######
  d$rowid <- seq.int(nrow(d))
  d$Group <- droplevels(factor(d$Group, levels = as.character(chroms)))
  
  dinfo <- d[, colnames(d) %in% c("rowid", "Color", "pval", "Location", "Shape", "is_highlight", "is_annotate", "external_gene_name"), drop=FALSE]
  
  ###### Create position index ######
  subd <- d[, c("Variable", "Group", "pvalue", "rowid")]
  d_order <- subd[order(subd$Group, subd$Variable),]
  d_order$pos_index <- seq.int(nrow(d_order))
  
  ### Set up dataframe with position info ###
  subd <- d_order[, colnames(d_order)!="rowid"]
  maxRows <- by(subd, subd$Group, function(x) x[which.max(x$pos_index),])
  minRows <- by(subd, subd$Group, function(x) x[which.min(x$pos_index),])
  milimits <- do.call(rbind, minRows)
  malimits <- do.call(rbind, maxRows)
  lims <- merge(milimits, malimits, by="Group")
  names(lims) <- c("Color", "Varx", "px", "posmin", "Vary", "py", "posmax")
  lims$av <- (lims$posmin + lims$posmax)/2
  lims$shademap <- rep(c("shade_ffffff","shade_ebebeb"), each=1, length.out=nrow(lims))
  
  ###### Color palettes ######
  ### Fill
  nvarcolors <- nlevels(factor(lims$Color))
  base_color <- c(rep(x=c(color1, color2), 
                      length.out=nvarcolors, each=1), 
                  "#FFFFFF", "#EBEBEB")
  names(base_color) <- c(levels(factor(lims$Color)), 
                         "shade_ffffff", "shade_ebebeb")
  ### Color
  if("Color" %in% names(d)){
    ### Color by Color column
    if(!missing(groupcolors)){
      dcols <- c(rep(x=c(color1, color2), length.out=nvarcolors, each=1), "#FFFFFF", "#EBEBEB")
      names(dcols) <-c(levels(factor(lims$Color)), "shade_ffffff", "shade_ebebeb")
      topcols <- c(dcols, groupcolors)
      bottomcols <- c(dcols, groupcolors)
    } else {
      
      ### Top Colors
      ngroupcolors <- nlevels(factor(d$Color[d$Location=="Top"]))
      if(ngroupcolors > 15){
        topcols <- Turbo(out.colors=ngroupcolors)
      } else {
        pal <- c("#009292", "#920000", "#490092", "#db6d00", "#24ff24", 
                 "#ffff6d", "#000000", "#006ddb", "#004949","#924900", 
                 "#ff6db6", "#6db6ff","#b66dff", "#ffb6db","#b6dbff")
        topcols <- pal[1:ngroupcolors]
      }
      names(topcols) <- levels(factor(d$Color[d$Location=="Top"]))
      
      ### Bottom Colors
      ngroupcolors <- nlevels(factor(d$Color[d$Location=="Bottom"]))
      if(ngroupcolors > 15){
        if (!requireNamespace(c("RColorBrewer"), quietly = TRUE)==TRUE) {
          stop("Please install RColorBrewer to add color attribute for more than 15 colors.", call. = FALSE)
        } else {
          bottomcols <- Turbo(out.colors=ngroupcolors)
        }
      } else {
        pal <- c("#009292", "#920000", "#490092", "#db6d00", "#24ff24", 
                 "#ffff6d", "#000000", "#006ddb", "#004949","#924900", 
                 "#ff6db6", "#6db6ff","#b66dff", "#ffb6db","#b6dbff")
        bottomcols <- pal[1:ngroupcolors]
      }
      names(bottomcols) <- levels(factor(d$Color[d$Location=="Bottom"]))
    }
  } else {
    ### Color by 'Group' instead
    names(d_order)[names(d_order)=="Group"] <- "Color"
    topcols <-c(rep(x=c(color1, color2), length.out=nvarcolors, each=1), "#FFFFFF", "#EBEBEB")
    names(topcols) <-c(levels(factor(lims$Color)), "shade_ffffff", "shade_ebebeb")
    bottomcols <-c(rep(x=c(color1, color2), length.out=nvarcolors, each=1), "#FFFFFF", "#EBEBEB")
    names(bottomcols) <-c(levels(factor(lims$Color)), "shade_ffffff", "shade_ebebeb")
  }
  
  ###### Theme options ######
  backpanel1 <- ifelse(background=="white", "NULL", 
                       "geom_rect(data = lims, aes(xmin = posmin-.5, xmax = posmax+.5, ymin = yaxismin1, ymax = Inf, fill=factor(shademap)), alpha = 0.5)" )
  backpanel2 <- ifelse(background=="white", "NULL", 
                       "geom_rect(data = lims, aes(xmin = posmin-.5, xmax = posmax+.5, ymin = yaxismin2, ymax = Inf, fill=factor(shademap)), alpha = 0.5)" )
  
  ###### Start plotting ######
  d_order <- merge(d_order, dinfo, by="rowid")
  
  ####### TOP PLOT ######
  p1 <- ggplot() + 
    eval(parse(text=backpanel1))
  
  ### Add shape info if available
  # if("Shape" %in% names(d)){
  #   p1 <- p1 + 
  #     geom_point(data=d_order[d_order$Location=="Top",], 
  #                aes(x=pos_index, 
  #                    y=pval, 
  #                    color=Color, 
  #                    shape=factor(Shape)), 
  #                alpha=opacity,
  #                size=0.5)
  # } else {
  #   p1 <- p1 + 
  #     geom_point(data=d_order[d_order$Location=="Top",], 
  #                aes(x=pos_index, 
  #                    y=pval, 
  #                    color=Color), 
  #                alpha=opacity,
  #                size=0.5)
  # }
  # if no shape info is added.
  p1 <- p1 + 
    geom_point(data=d_order[d_order$Location=="Top",],
               aes(x=pos_index,
                   y=pval,
                   color=Color),
               alpha=opacity,
               size=0.2)
  
  p1 <- p1 + 
    scale_x_continuous(breaks=lims$av, 
                       labels=lims$Color, 
                       expand=c(0,0))
  if(grpblocks==TRUE){
    if(freey==TRUE){
      print("Sorry, drawing grpblocks with freey=TRUE is currently unsupported and will be ignored.")
    } else {
      p1 <- p1 + 
        geom_rect(data = lims, 
                  aes(xmin = posmin-.5, 
                      xmax = posmax+.5, 
                      ymin = -Inf, 
                      ymax = min(d_order$pval), 
                      fill=as.factor(Color)), 
                  alpha = 1)
    }
  }
  
  p1 <- p1 +
    theme(panel.grid.minor.x = element_blank(), 
          panel.grid.major.x=element_blank(), 
          axis.title.x=element_blank(), 
          legend.position="top", 
          legend.title=element_blank())
  
  
  ###### BOTTOM PLOT ######
  p2 <- ggplot() + 
    eval(parse(text=backpanel2))
  
  ### Add shape info if available
  # if("Shape" %in% bottomn){
  #   p2 <- p2 + 
  #     geom_point(data=d_order[d_order$Location=="Bottom",], 
  #                aes(x=pos_index, 
  #                    y=pval, 
  #                    color=Color, 
  #                    shape=factor(Shape)), 
  #                alpha=opacity,
  #                size=0.2)
  # } else {
  #   p2 <- p2 + 
  #     geom_point(data=d_order[d_order$Location=="Bottom",], 
  #                aes(x=pos_index, 
  #                    y=pval, 
  #                    color=Color), 
  #                alpha=opacity,
  #                size=0.2)
  # }
  # add in color coding
  p2 <- p2 + 
    geom_point(data=d_order[d_order$Location=="Bottom",], 
               aes(x=pos_index, 
                   y=pval, 
                   color=Color), 
               alpha=opacity,
               size=0.2)
  p2 <- p2 + 
    scale_x_continuous(breaks=lims$av, 
                       labels=lims$Color, 
                       expand=c(0,0))
  if(grpblocks==TRUE){
    p2 <- p2 + 
      geom_rect(data = lims, 
                aes(xmin = posmin-.5,
                    xmax = posmax+.5, 
                    ymin = -Inf, 
                    ymax = min(d_order$pval), 
                    fill=as.factor(Color)), 
                alpha = 1)
  }
  
  p2 <- p2 + 
    theme(panel.grid.minor.x = element_blank(), 
          panel.grid.major.x=element_blank(), 
          axis.title.x=element_blank(), 
          legend.position="bottom", 
          legend.title=element_blank())
  
  
  ###### Legends ######
  p1 <- p1 + 
    scale_colour_manual(name = "Color", 
                        values = topcols) + 
    scale_fill_manual(values = base_color) +
    guides(fill="none")
  
  p2 <- p2 +
    scale_colour_manual(name = "Color", 
                        values = bottomcols) + 
    scale_fill_manual(values=base_color) +
    guides(fill="none")
  
  if(!("Color" %in% names(d))){
    p1 <- p1 + guides(color="none", fill="none")
    p2 <- p2 + guides(color="none", fill="none")
  }
  
  ####### Highlighting ######
  if(!missing(highlight_var)){
    if("Shape" %in% topn){
      p1 <- p1 + 
        geom_point(data=d_order[d_order$Variable %in% highlight_var & d_order$Location=="Top", ], 
                   aes(x=pos_index, 
                       y=pval),
                   size = 0.6,
                   #shape=Shape), 
                   colour=highlighter)
      p1 <- p1 + 
        guides(shape = guide_legend(override.aes = list(colour = "black")))
    } else {
      p1 <- p1 + 
        geom_point(data=d_order[d_order$Variable %in% highlight_var & d_order$Location=="Top", ], 
                   aes(x=pos_index, 
                       y=pval),
                   size =0.6,
                   colour=highlighter)
    }
    if("Shape" %in% bottomn){
      p2 <- p2 + 
        geom_point(data=d_order[d_order$Variable %in% highlight_var & d_order$Location=="Bottom", ], 
                   aes(x=pos_index, 
                       y=pval),#                       shape=Shape), 
                   size =0.6,
                   colour=highlighter)
      p2 <- p2 + guides(shape = guide_legend(override.aes = list(colour = "black")))
    } else {
      p2 <- p2 + 
        geom_point(data=d_order[d_order$Variable %in% highlight_var & d_order$Location=="Bottom", ], 
                   aes(x=pos_index, 
                       y=pval),
                   size =0.6,
                   colour=highlighter)
    }
  }
  if(!missing(highlight_p)){
    if("Shape" %in% topn){
      p1 <- p1 + 
        geom_point(data=d_order[d_order$pvalue < highlight_p[1] & d_order$Location=="Top", ], 
                   aes(x=pos_index, 
                       y=pval), #                       shape=Shape), 
                   colour=highlighter)
      p1 <- p1 + guides(shape = guide_legend(override.aes = list(colour = "black")))
    } else {
      p1 <- p1 + 
        geom_point(data=d_order[d_order$pvalue < highlight_p[1] & d_order$Location=="Top", ], 
                   aes(x=pos_index, 
                       y=pval), 
                   colour=highlighter)
    }
    if("Shape" %in% bottomn){
      p2 <- p2 + 
        geom_point(data=d_order[d_order$pvalue < highlight_p[2] & d_order$Location=="Bottom", ], 
                   aes(x=pos_index, 
                       y=pval), #                     shape=Shape), 
                   colour=highlighter)
      p2 <- p2 + guides(shape = guide_legend(override.aes = list(colour = "black")))
    } else {
      p2 <- p2 + geom_point(data=d_order[d_order$pvalue < highlight_p[2] & d_order$Location=="Bottom", ], 
                            aes(x=pos_index, 
                                y=pval), 
                            colour=highlighter)
    }
  }
  
  ####### Annotations ######
  if(!missing(annotate_p)){
    if (!requireNamespace(c("ggrepel"), quietly = TRUE)==TRUE) {
      print("Consider installing 'ggrepel' for improved text annotation")
      p1 <- p1 + geom_text(data=d_order[d_order$pvalue < annotate_p[1] & d_order$Location=="Top",], 
                           aes(pos_index,
                               pval,
                               label=Variable))
      
      p2 <- p2 + geom_text(data=d_order[d_order$pvalue < annotate_p[2] & d_order$Location=="Bottom",], 
                           aes(pos_index,
                               pval,
                               label=Variable))
    } else {
      p1 <- p1 + ggrepel::geom_text_repel(data=d_order[d_order$pvalue < annotate_p[1] & d_order$Location=="Top",], 
                                          aes(pos_index,
                                              pval,
                                              label=Variable))
      p2 <- p2 + ggrepel::geom_text_repel(data=d_order[d_order$pvalue < annotate_p[2] & d_order$Location=="Bottom",], 
                                          aes(pos_index,
                                              pval,
                                              label=Variable))
    }
  }
  
  # add gene labels 
  if(!missing(annotate_var)){
    if (!requireNamespace(c("ggrepel"), quietly = TRUE)==TRUE){
      print("Consider installing 'ggrepel' for improved text annotation")
      p1 <- p1 + geom_label_repel(data=d_order[d_order$Variable %in% annotate_var & d_order$Location=="Top",], 
                                  aes(pos_index,
                                      pval,
                                      label=external_gene_name),
                                  max.overlaps = Inf,
                                  seed = 42,
                                  box.padding = 0.5, 
                                  size = 3,
                                  #direction="y",
                                  ylim = c(4,8))
      p2 <- p2 + geom_label_repel(data=d_order[d_order$Variable %in% annotate_var & d_order$Location=="Bottom",], 
                                  aes(pos_index,
                                      pval,
                                      label=external_gene_name),
                                  max.overlaps = Inf,
                                  seed = 42,
                                  box.padding = 0.5,
                                  size = 3,
                                  #direction="y",
                                  ylim = c(-4,-8))
    } else {
      p1 <- p1 + ggrepel::geom_label_repel(data=d_order[d_order$Variable %in% annotate_var & d_order$Location=="Top",], 
                                           aes(pos_index,
                                               pval,
                                               label=external_gene_name),
                                           max.overlaps = Inf,
                                           seed = 42,
                                           # direction="y",
                                           box.padding = 0.5,
                                           size = 3,
                                           ylim = c(2,8))
      
      p2 <- p2 + ggrepel::geom_label_repel(data=d_order[d_order$Variable %in% annotate_var & d_order$Location=="Bottom",], 
                                           aes(pos_index,
                                               pval,
                                               label=external_gene_name),
                                           max.overlaps = Inf,
                                           # direction="y",
                                           seed = 42,
                                           box.padding = 0.5,
                                           size = 3,
                                           ylim = c(-2,-8))
    }
  }
  
  ###### Additional attributes ######
  # Add pvalue threshold line
  if(!missing(tline)){
    for(i in 1:length(tline)){
      p1 <- p1 + geom_hline(yintercept = tredline[i], colour="red", linetype = "dashed", alpha = 0.5)
    }
  }
  if(!missing(bline)){
    for(i in 1:length(bline)){
      p2 <- p2 + geom_hline(yintercept = bredline[i], colour="red", linetype = "dashed", alpha = 0.5)
    }
  }
  
  # Add title and y axis title
  p1 <- p1 + ylab(yaxislab1)
  p2 <- p2 + ylab(yaxislab2)
  
  # Formatting
  if(grpblocks==TRUE){
    p1 <- p1+ theme(text = element_text(size=8),
                    axis.text.x = element_text(vjust=5.5, size = 6),
                    axis.ticks.x = element_blank()) +
      ylim(c(yaxismin1, yaxismax1))
    p2 <- p2 + scale_y_reverse(limits=c(yaxismax2, yaxismin2)) + 
      theme(text = element_text(size=8),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank())
  } else {
    p1 <- p1 + theme(text = element_text(size = 8),
                     axis.text.x = element_text(vjust=5.5, size =6),
                     axis.ticks.x = element_blank()) + 
      scale_y_continuous(limits=c(0, 8))
    p2 <- p2 + scale_y_reverse(limits=c(8,0)) + 
      theme(text = element_text(size = 8),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank())
  }
  
  
  if(background=="white"){
    p1 <- p1 + theme(panel.background = element_rect(fill="white"))
    p2 <- p2 + theme(panel.background = element_rect(fill="white"))
  }
  
  # position the two plots
  p1 <- p1 + theme(plot.margin=unit(c(5.5,5.5,-3,5.5), unit = "pt"))
  p2 <- p2 + theme(plot.margin=unit(c(-3,5.5,5.5,5.5), unit = "pt"))
  
  if(rotatelabels==TRUE){p1 <- p1 + theme(axis.text.x = element_text(angle=labelangle))}
  
  
  # g_legend<-function(a.gplot){
  #   tmp <- ggplot_gtable(ggplot_build(a.gplot))
  #   leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  #   legend <- tmp$grobs[[leg]]
  #   return(legend)}
  # 
  # mylegend <- g_legend(p1 + theme(legend.position='top'))
  # 
  p <- grid.arrange(
    arrangeGrob(
      #mylegend,
      p1 + theme(legend.position="none"),
      p2 + theme(legend.position="none"),
      padding = 0),
    top = textGrob(toptitle, gp=gpar(fontsize=12,font=1)),
    heights=c(hgtratio, hgtratio))
  
  ###### Save ######
  # print(paste0("Saving plot to ", file, ".", type))
  # p <- grid.arrange(arrangeGrob(p1, top=toptitle),
  #                   arrangeGrob(p2, bottom=bottomtitle),
  #                   padding=0,
  #                   heights=c(hgtratio, 1-hgtratio))
  # ggsave(p, filename=paste0(file, ".", type), dpi=res, units="in", height=hgt, width=wi)
  return(p)
}


## function to apply adapted_emirror to all datasets
manhattan_continuous <- function(phenotype, biotype = "promoters", category = "continuous"){

  for (region in all_regions){
    
    celltypespresent <- readRDS(paste0(baseDir,"/cleaned_mvals/cleaned_mvals_episcore/episcore_celltypes_present/",region, "celltypes_present.rds"))
    celltypespresent <- c("bulk", celltypespresent)
    
    for (celltype in celltypespresent) {
      
      # load df with all protein coding results 
      file_results <- paste0(baseDir, "/dm_results/", category, "/raw/protein_coding/raw_protein_coding", phenotype, "_",  celltype, "_",biotype, "_", region, ".csv")
      if(file.exists(file_results)){
        all_results <- read.csv(file_results)
        all_results$ID <- all_results$ensembl_gene_id
      } 
      # load df with significant protein coding results 
      file_results <- paste0(baseDir, "/dm_results/", category, "/sig/protein_coding/protein_coding_sig_", phenotype, "_",  celltype, "_",biotype, "_", region, ".csv")
      if(file.exists(file_results)){
        sig_res <- read.csv(file_results)
        sig_res$ID <- sig_res$ensembl_gene_id
      # get top 20 protein coding highest and lowest logFC results to highlight in manhattan plot
      pos_log_fc <- sig_res[sig_res$logFC > 0,]
      pos_log_fc <- pos_log_fc[order((pos_log_fc$logFC), decreasing = TRUE),]
      pos_log_fc <- pos_log_fc[1:20,]
      neg_log_fc <- sig_res[sig_res$logFC < 0,]
      neg_log_fc <- neg_log_fc[order(abs(neg_log_fc$logFC), decreasing = TRUE),]
      neg_log_fc <- neg_log_fc[1:20,]
      annotate_top20_log_fc <- rbind(neg_log_fc, pos_log_fc)
      
      # mark significant protein coding results and for top20 significant protein coding results add annotation notice
      all_results <- all_results %>% mutate(is_annotate=ifelse(ensembl_gene_id %in% annotate_top20_log_fc$ensembl_gene_id, "yes", "no"))
      # get promoter annotation
      df2 <- anno_promoters[,c("Start", "Chromosome", "ID")]
      
      # merge all results with promoter annotation to add in Chromosme information, start and End
      plot_df <- merge(all_results, df2, by = "ID", all.x=TRUE)
      
      # make sure that chromosome start site is numeric
      plot_df$Start = as.numeric(as.vector(plot_df$Start))
      plot_df = plot_df[order(plot_df$Chromosome,plot_df$Start),]
      
      # important for plotting: add in a cumulative score ens_cum for the position of each Ens-ID
      plot_df <- plot_df %>%
        # Compute chromosome size
        group_by(Chromosome) %>% 
        dplyr::summarise(chr_len=max(Start)) %>% 
        
        # Calculate cumulative position of each chromosome
        mutate(tot=cumsum(chr_len)-chr_len) %>%
        dplyr::select(-chr_len) %>%
        
        # Add this info to the initial dataset
        left_join(plot_df, ., by=c("Chromosome"="Chromosome")) %>%
        
        # Add a cumulative position of each ens id 
        arrange(Chromosome, Start) %>%
        mutate(ens_cum = Start+tot)
      
      plot_df$Chromosome <- str_replace(plot_df$Chromosome, "chr", "")
      
      positive <- plot_df  %>% filter(plot_df$logFC > 0) %>% mutate(group = "positive")
      negative <- plot_df  %>% filter(plot_df$logFC < 0) %>% mutate(group = "negative") 
      
      pos <- positive %>% dplyr::select(c("ens_cum", "adj.P.Val", "Chromosome", "is_annotate", "external_gene_name"))
      colnames(pos) <- c("Variable", "pvalue", "Group", "is_annotate", "external_gene_name")
      
      neg <- negative %>% dplyr::select(c("ens_cum", "adj.P.Val", "Chromosome", "is_annotate", "external_gene_name"))
      colnames(neg) <- c("Variable", "pvalue", "Group", "is_annotate", "external_gene_name")
      
      
      #-------------------------------------------------------------------------------
      # mirror manhattan plot with labels and colored top20 protein coding results
      #-------------------------------------------------------------------------------
      annotate <- plot_df[plot_df$is_annotate == "yes",]$ens_cum
      p <- adapted_emirror(top=pos, bottom=neg, tline = 0.05, bline = 0.05, color1="grey40", highlighter = "#07D4FF",
                           color2="grey80", rotatelabels = FALSE,labelangle = 45, background="white",opacity= 0.5, hgtratio=4, highlight_var = annotate, annotate_var = annotate,
                           toptitle =  "")

      ggsave(
        paste0(saveDir, "plots/final_figures/Fig4_manhattan_mirror_", category, "_", phenotype, "_",region, "_", celltype, "_", biotype, ".png"),
        p,
        width = 6.85,
        height = 6.85,
        dpi = 200,
        units = 'in'
      )
      
      pdf(file = paste0(saveDir, "plots/final_figures/Fig4_manhattan_mirror_", category, "_", phenotype, "_",region, "_", celltype, "_", biotype, ".pdf"),
      width = 6.85,
      height = 6.85)
      print(p)
      dev.off()
        
      }
      }
  }
}




## run manhattan mirror plot
manhattan_continuous(phenotype = "niaaaascore", category= "continuous")
manhattan_continuous(phenotype = "niaaabscore", category= "continuous")
manhattan_continuous(phenotype = "int_adseverityscore", category= "continuous")
manhattan_continuous(phenotype = "braak_for_lb_orig", category= "continuous")
manhattan_continuous(phenotype = "mvl_score", category= "continuous")
manhattan_continuous(phenotype = "tdp43", category= "continuous")
