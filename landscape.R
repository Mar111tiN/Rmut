library(tidyverse)
library(patchwork)

############# WORKFLOW FOR LANDSCAPE DATA ################
###### consists of these data munging steps:
###### -  read.mut(file):  load edited filter1 list from somVar pipeline
###### -  filter2(): apply filter2 using a filter setting provided as excel file
###### -  remove.dups: remove duplicate muts per gene (at different positions)
###### -  count.gene(): 
###### -  format4landscape(): condense the filter info for landscape function
###### -  plot.landscape(): run landscape vis in ggplot2


#############################################################
##### remove duplicate entries ##############################
remove.dups <- function(df) {
  ## remove duplicate entries - same sample gene - different position
  ## add genesPerSample
  condensed.df <- df %>% 
    # select important columns
    select(sample, Seq, Chr, Start, Alt, Gene, Func, isCandidate, AMLDriver, ChipFreq, TVAF, NVAF) %>%
    # remove duplicate entries (happened for KRAPK5-10, weirdly)
    distinct(sample, Seq, Chr, Gene, Start, Alt, .keep_all = T) %>% 
    # for each patient, only keep the last sequenced mutations
    group_by(sample) %>%
    filter(rank(desc(Seq)) == min(rank(desc(Seq)))) %>% 
    # group for condensing the duplicate mutations
    group_by(sample, Seq, Gene) %>% 
    # add duplication count (measure of dirtyness?)
    mutate(multipleCount = n()) %>% 
    # per group, only recover the ones with highest TVAF and existing ChipFreq
    filter(rank(desc(ChipFreq)) == min(rank(desc(ChipFreq)))) %>%  
    filter(rank(desc(TVAF)) == min(rank(desc(TVAF)))) %>% 
    ungroup() %>% 
    group_by(sample) %>% 
    mutate(genesPerSample = n()) %>% 
    ungroup()
  print(glue('[remove.dups] Dropped {nrow(df) - nrow(condensed.df)} of {nrow(df)} mutations --> {nrow(condensed.df)}'))
  return(condensed.df)
}

####################################################
##### format4landscape #############################
count.gene <- function(df, min.count=3, min.count.cand=1, min.count.driver=2, hotspot.count=2, min.count.hotspot=1) {
  count.df <- df %>% 
    remove.dups() %>% 
    select(sample, Gene, Func, isCandidate, AMLDriver, TVAF, ChipFreq) %>% 
    # translate Func into consequence for color scheme
    # add sampleCount
    group_by(sample) %>% 
    mutate(sampleCount = n()) %>% 
    # add geneCount
    group_by(Gene) %>% 
    mutate(geneCount = n()) %>% 
    # filter by minimum gene count or candidate
    filter(geneCount >= min.count | (isCandidate == 1  & geneCount >= min.count.cand) | (AMLDriver == 1 & geneCount >= min.count.driver) | ChipFreq >= hotspot.count  & geneCount >= min.count.hotspot) %>% 
    ungroup() %>% 
    # reorder categorical Gene depending on gene_count 
    mutate(geneOrdered = fct_reorder(as.factor(Gene), geneCount)) %>%
    # create gene_rank as numerical value and convert to summable exponent
    mutate(geneRank = as.numeric(geneOrdered)) %>% 
    group_by(sample) %>% 
    # add sampleCount after reduction of minimum counts
    mutate(sampleCountCut = n()) %>% 
    ungroup()
  # output
  total.count <- nrow(df)
  count.count <- nrow(count.df)
  drop.count <- total.count - count.count
  print(glue('[count.gene (min.count = {min.count})] Dropped {drop.count} of {total.count} mutations --> {count.count}'))
  return(count.df)
}

############### HELPER for LANDSCAPE PLOT ##################
arrange4landscape <- function(df.counted) {
  return(
    df.counted %>% 
      # add helper geneBinary for summing up the genes in ranking order
      mutate(geneBinary = 2^geneRank) %>% 
      # group by sample and rank patients according to ranked gene sum 
      group_by(sample) %>% 
      # summing up the binaryGenes for correct sorting of the patients
      mutate(geneScore = sum(geneBinary)) %>% 
      group_by(Func) %>% 
      mutate(FuncCount = n()) %>% 
      ungroup() %>% 
      # create an ordered sample factor to be used by the tile geom
      mutate(sampleOrdered = fct_reorder(as.factor(sample), -geneScore)) %>% 
      select(sampleOrdered, geneOrdered, Func, isCandidate, AMLDriver, ChipFreq, geneRank, sampleCount, sampleCountCut)
  )
}

#######################################################################
################# LANDSCAPE PLOT FUNCTION #############################
plot.waterfall <- function(df.arranged, color.code=c( ### SET THE COLOURS FOR THE LANDSCAPE GRAPH
  `Missense Mutation` = "orangered3",
  `Splicing` = "deepskyblue4",
  `Nonsense Mutation` = "grey57",
  `Frameshift Indel` = "firebrick4",
  `Inframe Indel` = "darkgoldenrod1",
  `Multiple` = "yellow",
  `not Exonic` = "gray"
)) {
  return(df.arranged %>% 
           ggplot(aes(
             x = sampleOrdered,
             y = geneOrdered,
             fill = Func
           ))+
           geom_tile() +
           labs(
             x = NULL,
             y = "Genes",
             fill = "Mutation"
           ) +
           scale_fill_manual(values = color.code) +
           theme(
             plot.margin = margin(0,0,0,0),
             panel.grid = element_blank(),
             axis.ticks = element_blank(),
             axis.text.x = element_text(angle=90, vjust=0.5),
             # axis.text.x = element_blank(),
             axis.text.y = element_text(
               hjust = 0.5
             ),
             axis.title.y = element_blank()
           )
  )
}
################# sampleHisto PLOT FUNCTION (horizontal plot) #############################
plot.sampleHisto <- function(df.arranged) {
  return(df.arranged %>% 
           ggplot() +
           geom_col(aes(sampleOrdered, sampleCount)) +
           geom_col(aes(sampleOrdered, sampleCountCut), fill="red") +
           scale_y_reverse(
             name = quote(frac(mutations, sample)),
             expand= c(0,0),
             breaks = c(20,100,200,500,1000,2000)
           ) +
           theme(
             plot.margin = margin(0,0,0,0),
             axis.text.x = element_blank(),
             axis.ticks.x = element_blank(),
             axis.title.x = element_blank()
           )
  )
}
################# geneHisto PLOT FUNCTION (vertical plot) #############################
plot.geneHisto <- function(df.counted) {
  return(
    df.counted %>% 
      select(Gene, geneRank, geneCount) %>% 
      group_by(Gene) %>%  
      summarise(geneCount = unique(geneCount), geneRank = first(geneRank)) %>% 
      arrange(-geneRank) %>% 
      ggplot() +
      geom_col(aes(fct_reorder(Gene, geneCount), geneCount)) +
      labs(x = "Gene count") +
      theme(
        plot.margin = margin(0,0,0,0),
        axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()
      ) +
      scale_y_reverse(
        expand = c(0,0),
        breaks = c(0,5,10,15)
      ) +
      coord_flip()
  )
}

######## MAKE COMBINED LANDSCAPE PLOT
plot.landscape <- function(df, ...) {
  # create subdata
  df.counted <- df %>% 
    count.gene(...) 
  df.arranged <- df.counted %>% 
    arrange4landscape()
  # create subplots
  waterfall.plot <- plot.waterfall(df.arranged)
  sampleHisto.plot <- plot.sampleHisto(df.arranged)
  geneHisto.plot <- plot.geneHisto(df.counted)
  return(
    geneHisto.plot + waterfall.plot + plot_spacer() + sampleHisto.plot + 
      plot_layout(ncol = 2, widths = c(1,10), heights = c(10,1))
  )
}