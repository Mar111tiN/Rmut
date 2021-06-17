# install.packages('tidyverse')
library(tidyverse)
library(readxl)

# specialized parse function
read.mut <- function(file.path, factor.cols = c('Seq', "Daniel-ID", 'Diagnosis', 'Chr', 'Func')) {
  parser.cols <- cols_condense(spec_tsv(file.path, guess_max = 10000))
  # change parse options for certain cols
  for (col in factor.cols) { parser.cols$cols[[col]] = col_factor() }
  # parse file into filter1 (tibble)
  filter1 <- read_tsv(file.path, col_types = parser.cols)
  return(filter1)
}


condense.func <- function(df) {
  df <- df %>% 
    mutate(Func = str_replace(Func, "upstream|downstream|intergenic|^synSNV|;synSNV|ncRNA_splicing|ncRNA|unknown|UTR3|UTR5", ""))%>% 
    mutate(Func = str_replace(Func, "intronic|exonic", "")) %>% 
    mutate(Func = str_replace(Func, "_", "")) %>% 
    mutate(Func = str_replace(Func, "^;+", "")) %>% 
    mutate(Func = str_replace(Func, ";+$", "")) %>% 
    mutate(Func = str_replace(Func, ";;", ";"))%>% 
    mutate(Func = str_replace(Func, "upstream|downstream|intergenic|^synSNV|;synSNV|ncRNA_splicing|ncRNA|unknown|UTR3|UTR5", ""))%>% 
    mutate(Func = str_replace(Func, "intronic|exonic", "")) %>% 
    mutate(Func = str_replace(Func, "_", "")) %>% 
    mutate(Func = str_replace(Func, "^;+", "")) %>% 
    mutate(Func = str_replace(Func, ";+$", "")) %>% 
    mutate(Func = str_replace(Func, ";;", ";"))

    # chance columns
  trans.consequence <- c(
    `not Exonic`= "not Exonic",
    Multiple = "Multiple",
    `nonsynSNV` = 'Missense Mutation',
    splicing = 'Splicing',
    stoploss = 'Nonsense Mutation',
    `frameshift deletion` = 'Frameshift Indel',
    `nonframeshift deletion` = "Inframe Indel",
    startloss = 'Nonsense Mutation',
    stopgain = 'Nonsense Mutation',
    `frameshift insertion` = 'Frameshift Indel',
    `nonframeshift insertion` = "Inframe Indel"
  )

  df <- df %>%
    mutate(Func = case_when(Func == "" ~ "not Exonic", str_detect(Func, ";") ~ "Multiple", TRUE ~ Func)) %>%
    mutate(Func = trans.consequence[Func])

  return(df)
}


###############################################################
######### FILTER2 ##############################################

# helper function to load the desired filter as a tibble of dim c(1,N)
get.filter <- function(filter.file, sheet = 'AMLMono7', stringency='moderate', ...) {
  # load excel file and convert 
  settings <- suppressMessages(read_excel(filter.file, sheet = sheet, n_max = 4) %>%
                                 mutate_at(vars(variantT, Tdepth, PONAltNonZeros, HDRcount), list(as.integer)) %>% 
                                 mutate_at(vars(EBscore, PONRatio, FisherScore, TVAF, TVAF4Cand, NVAF, VAFSim, ClinScore), list(as.numeric)))
  rows <- 1:4
  if (!grepl("filter", stringency)) stringency <- paste("filter2-", stringency, sep='')
  names(rows) <- c("filter1", "filter2-loose", "filter2-moderate", "filter2-strict")
  # restrict to the filter rows (other rows are only for users)
  f = settings[rows[stringency], ,drop=TRUE]
  # update the filter by additional params
  mod.params <- as.vector(list(...))
  for (param in names(mod.params)) {
    f[param] <- mod.params[param]
  }
  return(f)
}


### filter function #################
# uses the filter tibble f and runs through the filter steps
# the ClinScore rescue is performed via a full outer join of the filtered muts and the ClinScore+ muts
filter2 <- function(df, filter.file='', sheet='AMLMono7', stringency='moderate', use.validate=F, drop.ambiguous=F, skip.samples=c(), ...) {
  ### get filter settings with overloaded additional params
  f <- get.filter(filter.file, sheet=sheet, stringency=stringency, ...)

  # drop nonExonic muts
  df <- df  %>% 
    filter(Func != "not Exonic")
    
  ## VALIDATION FILTER #############
  if (use.validate) {
    if (drop.ambiguous) {
      df.validated <- df %>% 
        filter(!Call %in% c('G', 'F', 'A'))
    } else {
      df.validated <- df %>% 
        filter(!Call %in% c('G', 'F'))
    }
  } else {
    df.validated <- df
  }


  
  ## MAIN DATA PIPE ####
  filter2.true <- df.validated %>%
    # skip samples
    filter(!sample %in% skip.samples) %>% 
    # TumorDepth
    filter(TR2 > f$variantT & Tdepth > f$Tdepth)%>%
    # EB_PON
    filter((EBscore >= f$EBscore & PONRatio < f$PONRatio) | PONAltNonZeros < f$PONAltNonZeros) %>% 
    # TVAF
    filter(((isCandidate == 1 | AMLDriver == 1 | ChipFreq > 0) & TVAF >= f$TVAF4Cand) | TVAF >= f$TVAF) %>% 
    # NVAF
    filter(NVAF <= f$NVAF & TVAF >= f$VAFSim * NVAF) %>% 
    # Strandedness = Fisher & polarity
    filter(FisherScore <= f$FisherScore) %>% 
    filter((`TR2+` / TR2 <= f$strandPolarity) & (`TR2+` / TR2 >= 1 - f$strandPolarity)) %>% 
    filter(TumorHDRcount <= f$HDRcount & NormalHDRcount <= f$HDRcount)
  # rescued muts via ClinScore threshold
  filter2.rescue <- df %>% 
    filter(ClinScore >= f$ClinScore)
  filter2.merged <- suppressMessages(full_join(filter2.true, filter2.rescue))
  # output
  df.count <- nrow(df)
  merged.count <- nrow(filter2.merged)
  dropped <- df.count - merged.count
  resc.count <- merged.count - nrow(filter2.true)
  print(glue('[filter2-{stringency}] Dropped {dropped} of {df.count} mutations ({resc.count} rescued via ClinScore > {f$ClinScore}) --> {merged.count}'))
  return(filter2.merged)
}