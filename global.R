library(shiny)
library(gamlss)
library(magrittr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(Biostrings)
library(sangerseqR)

# functions copied from global.R
CreateSangs <- function(peakAmp, basecalls){
  # this function takes the peak amp matrix, and the input basecalls, and creates an 
  # unfiltered sanger sequencing dataframe (that I like to call sangs)
  sangs <- as.data.frame(peakAmp)
  names(sangs) <- c("A.area","C.area","G.area","T.area")
  
  # adding total area, and then calculating the percent area of each base
  sangs %<>% 
    mutate(Tot.area = A.area + C.area + G.area + T.area,
           A.perc = 100*A.area / Tot.area,
           C.perc = 100*C.area / Tot.area,
           G.perc = 100*G.area / Tot.area,
           T.perc = 100*T.area / Tot.area) 
  
  # adding on base calls
  sangs$base.call <- strsplit(x = toString(basecalls@primarySeq), split = "") %>%
    unlist
  
  # adding an index
  sangs$index <- seq_along(sangs$base.call)
  return(sangs)
}

GetGuideMatch <- function(guide, input.seq)  {
  # this function takes in the guide sequence, and the primary seq from a makeBaseCalls object
  # guide should be a DNAString object -- DNAString(guideinput)
  # input.seq = input.basecalls@primarySeq
  # returns guide.coord, a list with $match (forward / reverse), $start, and $end
  # and finds where the guide matches
  # it looks for both a forward and reverse match
  # in deciding if it's forward / reverse, it takes the match that is the longest
  # This function tries to extend the guide coordinates in case there is a partial match -- this might break
  # match to guide
  
  # alignments
  align.f  <- pairwiseAlignment(pattern=guide, subject=input.seq, type="overlap")
  
  guide.coord <- list(match = "forward", start = align.f@subject@range@start,
                      end = align.f@subject@range@start + align.f@subject@range@width-1)
  
  # check to see if a match was actually found
  if(align.f@pattern@range@width <= 1){
    stop("Failed to find a forward match of the gRNA -- should it be a reverse complement?")
    return(NA)
  }
  
  
  # check to see if the matched region == the guide length
  is.match.guide.length <- (guide.coord$end - guide.coord$start + 1) == length(guide)
  
  # otherwise figure out what part of the guide matched to where, and update the guide region to include the unmatched areas. 
  if(is.match.guide.length != TRUE){
    
    # find the diff from beginning
    startdiff <- 1 - align.f@pattern@range@start
    # find the diff from the end
    enddiff <- 20 - (align.f@pattern@range@start + align.f@pattern@range@width-1)
    # add these numbers to the guide.cord
    guide.coord$start <- guide.coord$start + startdiff
    guide.coord$end <- guide.coord$end + enddiff
    
    # check to make sure that I didn't do something bad
    # if(guide.coord$start < 1 | guide.coord$end > length(guide)){
    #   stop("Something went wrong in matching the guide to the sequence, partial match of the guide found.
    #          Error in extending the match to the full guide length.")
    # }
    
    
  }
  
  return(guide.coord)
}

GetNullDistModel <- function(sangs.filt, guide.coord) {
  # this function takes the sangs.filt dataframe, gathers values for the null distribution
  # for each base, and then fits a zero-adjusted gamma distribution to these values
  # it returns a list of data.frames, each with the parameters of a null distribution 
  # for each base
  
  sangs.filt %<>% filter(!(index %in% (guide.coord$start:guide.coord$end)) )
  
  nvals <- list()
  nvals$t <- sangs.filt %>% filter(base.call != "T") %>% dplyr::select(T.perc) %>% unlist()
  nvals$c <- sangs.filt %>% filter(base.call != "C") %>% dplyr::select(C.perc) %>% unlist()
  nvals$g <- sangs.filt %>% filter(base.call != "G") %>% dplyr::select(G.perc) %>% unlist()
  nvals$a <- sangs.filt %>% filter(base.call != "A") %>% dplyr::select(A.perc) %>% unlist()
  
  # Updated 3.26.19 to account for ultra clean sequencing
  replacement_zaga = c(rep(0, 989), 0.00998720389310502, 0.00998813447664401,0.009992887520785,
                       0.00999585366068316, 0.00999623914632598, 0.00999799013526835, 0.010001499423723,
                       0.0100030237039207, 0.0100045782875701, 0.0100048452355807, 0.0100049548867042)
  
  n.models <- lapply(nvals, FUN = function(x){
    set.seed(1)
    if((unique(x)[1] == 0 & length(unique(x)) == 1) |
       (unique(x)[1] == 0 & length(unique(x)) == 2 & table(x)[2] == 1))
    {x = replacement_zaga; message("Replacement vector used for low noise.")} # add noise if all 0s, or all 0s and one other value.
    tryCatch(gamlss((x)~1, family = ZAGA), error=function(e) # Progressively step up the mu.start if it fails
      tryCatch(gamlss((x)~1, family = ZAGA, mu.start = 1), error=function(e) 
        tryCatch(gamlss((x)~1, family = ZAGA, mu.start = 2), error=function(e) 
          tryCatch(gamlss((x)~1, family = ZAGA, mu.start = 3), error=function(e) # additional step added.
            gamlss((x)~1, family = ZAGA, mu.start = mean(x))
          )
        )
      )
    )
    # throws errors when a completely 0 vector
  })
  
  null.m.params <- lapply(n.models, FUN = function(x){
    mu <- exp(x$mu.coefficients[[1]])
    sigma <- exp(x$sigma.coefficients[[1]])
    nu.logit <- x$nu.coefficients[[1]]
    nu <- exp(nu.logit)/(1+exp(nu.logit))
    fillibens <-cor(as.data.frame(qqnorm(x$residuals, plot = FALSE)))[1,2]
    
    return(data.frame(mu= mu, sigma = sigma, nu = nu, fillibens = fillibens))
  })
  
  return(null.m.params)
}

CreateEditingDF <- function(guide.coord, guide, sangs, null.m.params){
  # this function creates the editing.df that contains the prob that the peak areas 
  # are not part of the noise distibution
  
  guide.df <- sangs[guide.coord$start:guide.coord$end,]
  
  # adding the guide sequence onto the data.frame
  guide.df$guide.seq <- guide %>% toString() %>% strsplit(. , "") %>% unlist
  
  
  # usage (params = null.m.params$t, perc = guide.df$T.perc)
  calcBaseProb <- function(params, perc){
    outprobs <- pZAGA(q = perc,
                      mu = params$mu,
                      sigma = params$sigma,
                      nu = params$nu,
                      lower.tail = FALSE)
  }
  
  guide.df$T.pval <- calcBaseProb(params = null.m.params$t, perc = guide.df$T.perc)
  guide.df$C.pval <- calcBaseProb(params = null.m.params$c, perc = guide.df$C.perc)
  guide.df$G.pval <- calcBaseProb(params = null.m.params$g, perc = guide.df$G.perc)
  guide.df$A.pval <- calcBaseProb(params = null.m.params$a, perc = guide.df$A.perc)
  
  guide.df$guide.position <- seq_along(guide.df$guide.seq)
  
  editing.df <- guide.df  
  return(editing.df)
}

GenSignalNoise <- function(sangs.filt) {
  catch <- sangs.filt %>% dplyr::select(A.area:T.area, A.perc:T.perc, base.call, index)
  catch %<>% gather( key = base, value = value,
                     A.area:T.area, A.perc:T.perc) %>%
    separate(col = base, into = c("base", "measure"))
  
  # splitting the catch dataframe into either signal and noise, and calculating the 
  # total noise area or total noise percent
  noise <- catch %>% 
    group_by(index, measure) %>% 
    filter(base != base.call) %>%
    summarize(noise = sum(value))
  signal <- catch %>%
    group_by(index, measure) %>%
    filter(base == base.call) %>% 
    summarize(signal = sum(value))
  
  signal.noise.df <- left_join(noise,signal) %>%
    gather(key = type, value = value, noise, signal)
}

CalculateTrim <- function(sangs, percent_noise, guide.coord) {
  # Custom trimming function
  # Idea: Use the percentage of the noise to the total signal to determine 
  # regions where the the reads get too noisy, where too noisy is defined by the
  # threshold, percent_noise 
  signal.noise.df <- GenSignalNoise(sangs)
  sangs.plot <- sangs %>% dplyr::select(index, base.call) 
  sangs.plot %<>% left_join(signal.noise.df)
  sangs.plot$type <- ordered(sangs.plot$type, levels = c("signal", "noise"))
  
  rolling_sang_percent <- sangs.plot %>% filter(measure == "perc", type == "noise") %>% dplyr::mutate(percent_3 = zoo::rollmean(value, k = 2, fill = NA)) %>% filter(percent_3 > percent_noise)
  
  peak_5 <- rolling_sang_percent %>% filter(index < guide.coord$start) %>% tail(n=1)
  peak_3 <- rolling_sang_percent %>% filter(index > guide.coord$end) %>% head(n=1)
  
  trim_left <- 20 # default start with 20 
  trim_right <- tail(sangs.plot, n=1)$index
  if (nrow(peak_5) > 0) {
    trim_left <- max(trim_left, peak_5$index)
  }
  if (nrow(peak_3) > 0) {
    trim_right <- peak_3$index
  }
  if ((trim_right - trim_left) < 100) {
    stop("Distance between autocalculated 5' trim and 3' trim location is less than 100 nucleotides. Read quality could be low, or alignmnet of the guide to the sequence might have been bad")
  }
  return(list(trim5 = trim_left, trim3 = trim_right))
}

CreateEditingSpread <- function(editing.df) {
  # Originally in server.R
  ### Reshape data
  edit.long <- editing.df %>% gather(key = focal.base, value = value, 
                                     A.area:T.area, A.perc:T.perc, T.pval:A.pval) %>%
    separate(col = focal.base, into = c("focal.base", "measure"))
  
  edit.spread <- edit.long %>% 
    spread(key = measure, value = value) 
}

CreateAvgBase <- function(sangs.filt) {
  #### Repeat code for getting avg.base from base.infoReactive
  # finding the average percent signal for each base
  sangs.filt %>% gather(key = focal.base, value = value, 
                        A.area:T.area, A.perc:T.perc) %>%
    separate(col = focal.base, into = c("focal.base", "measure")) %>% 
    spread(key = measure, value = value) %>% 
    filter(base.call == focal.base) %>% 
    group_by(focal.base) %>% 
    summarize(avg.percsignal = mean(perc),
              avg.areasignal = mean(area))
}

CalcEditR <- function(filename, guideseq, p.val.cutoff, default.trim, is.reverse) {
  # Returns editing information
  ## table - table of edits across the gRNA sequence
  ## figure - the editing efficiency plot on the "predicted editing" tab
  ## trim - the 5' and 3' trimming locations. If default.trim=TRUE, then the 5' is 20 and the 3' is 'Default' 
  
  input.seq = readsangerseq(filename)
  input.basecalls <- makeBaseCalls(input.seq)
  input.peakampmatrix <- peakAmpMatrix(input.basecalls)
  
  sangs <- CreateSangs(input.peakampmatrix, input.basecalls)
  
  # getting the sequence to match to
  filt.sequence <- sangs$base.call %>% paste(collapse = "") %>% DNAString()
  # this function finds where the guide matches
  guide <- DNAString(guideseq)
  if (is.reverse) {
    guide <- reverseComplement(guide)
  }
  guide.match <- GetGuideMatch(guide, filt.sequence)
  # Finding the index values
  guide.coord <- list(start = sangs[guide.match$start, "index"], end = sangs[guide.match$end, "index"])
  
  
  if (default.trim) {
    sangs.filt <- sangs %>% filter(index > 20)
    # removing crappy end
    peakTotAreaCutoff <- mean(sangs.filt$Tot.area)/10
    sangs.filt %<>% filter(Tot.area > peakTotAreaCutoff)
    trim.values <- list(trim5 = 20, trim3 = "Default")
  } else {
    trim.values <- CalculateTrim(sangs, 20, guide.coord)
    sangs.filt <- sangs[trim.values$trim5:trim.values$trim3, ]
  }
  
  # getting params for the different null models
  null.m.params <- GetNullDistModel(sangs.filt, guide.coord)
  editing.df <- CreateEditingDF(guide.coord, guide, sangs, null.m.params)
  edit.spread <- CreateEditingSpread(editing.df)
  avg.base <- CreateAvgBase(sangs.filt)
  
  # finding the model mu
  mul <- lapply(null.m.params, FUN = function(x){x$mu})
  mulvec <- c(a = mul$a, c = mul$c, g = mul$g, t = mul$t)
  
  color.cutoff = min(avg.base$avg.percsignal - mulvec)
  edit.color <- edit.spread %>% 
    mutate(adj.perc = {ifelse(perc >= color.cutoff,
                              100,
                              perc)
    } %>% as.numeric) %>%
    filter(pval < p.val.cutoff)
  
  edit.chart <- edit.spread %>%
    ggplot(aes(x = as.factor(index), y = focal.base)) + 
    geom_tile(data = edit.color, aes(fill = adj.perc)) + 
    geom_text(aes(label = round(perc, 0)), angle = 0, size = 5) +   
    guides(fill = FALSE) + 
    scale_fill_continuous(low = "#f7a8a8", high = "#9acdee") + 
    scale_x_discrete(position = "top", labels = editing.df$guide.seq) + 
    labs(x = NULL, y = NULL) + 
    theme(axis.ticks = element_blank(),
          axis.text=element_text(size=16),
          plot.title = element_text(hjust = 0, size = 16),
          plot.margin=unit(c(0,0,0,2), "cm"), #c(top, bottom, left, right)
          panel.background = element_rect(fill = "transparent",colour = NA), # or theme_blank()
          plot.background = element_rect(fill = "transparent",colour = NA)
    ) +
    coord_fixed(1)
  # count from the other direction if is reverse 
  if (is.reverse) {
    edit.spread$guide.position <- max(edit.spread$guide.position) + 1 - edit.spread$guide.position
  }
  
  return(list(table=edit.spread, figure=edit.chart, trim=trim.values))
}

CalcEditRBatch <- function(file.id, filenames, datapaths, guideseq, reverseYN) {
  p.val.cutoff <- 0.01
  default.trim <- FALSE
  is.reverse <- (reverseYN == "Y")
  
  file.matches <- grepl(paste0("^", file.id, "_"), filenames)
  if (sum(file.matches) != 1) {
    stop(paste("The provided file ID", file.id, "matched multiple files in the directory or could not be found"))
  }
  filename <- datapaths[file.matches]
  
  results <- CalcEditR(filename, guideseq, p.val.cutoff, default.trim, is.reverse)
  if (is.reverse) {
    results$table %<>% filter(guide.seq == "G") %>% filter(focal.base == "A")
  } else {
    results$table %<>% filter(guide.seq == "C") %>% filter(focal.base == "T")
  }
  results
}

