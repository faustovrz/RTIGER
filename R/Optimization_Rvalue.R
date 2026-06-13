#'
#' Find the otimum R value for a given data set
#'
#' @param object an RTIGER object
#' @param max_rigidity R values will be explored up the value given in this parameter. Default = 2^9
#' @param average_coverage For conservative results set it to the lowest average coverage of a sample in your experiment, or evne to the lowest average coverage in a (sufficiently large) region in one of your samples. The lower the value, the more conservative (higher) our estimates of the false positive segments rates. If it is not provided it will be computed as the average of all data points.
#' @param crossovers_per_megabase For conservative results set it to the highest ratio of a sample in your experiment. The higher the value, the more conservative (higher) our estimates of the false positive segments rates. If it is not provided it will be computed as the average of all samples.
#' @param save_it logical values if the results should be saved. Plots might be complicated to interpret. We suggest to read the manuscript to understand them (https://doi.org/10.1093/plphys/kiad191)
#' @param savedir if results are saved, in which directory.
#'
#' @return A value with the optimum rigidity for the data set.
#'
#' @usage optimize_R(object,
#' max_rigidity = 2^9, average_coverage = NULL, crossovers_per_megabase = NULL,
#' save_it = FALSE, savedir = NULL)
#'
#' @examples
#'
#' data("fittedExample")
#' bestR = optimize_R(myDat)
#'
#' @export optimize_R
#'

optimize_R = function(object,
                      max_rigidity = 2^9,
                      average_coverage = NULL,
                      crossovers_per_megabase = NULL,
                      seed = 1L,
                      n_obs = 1e4,
                      method = c("exact", "mc"),
                      save_it = FALSE,
                      savedir = NULL ){
  if(save_it & is.null(savedir)) stop("Please if you want to save the plots and results specify the path in savedir.\n")
  # method = "exact" (default): compute the FPR/FNR analytically by exponential-
  #   tilted FFT convolution of the per-marker log-likelihood-ratio increment.
  #   This is the manuscript's PAIRED statistic (Delta_{m,x,y} = logP(o,y)-logP(o,x)
  #   on the SAME observation o, Supp Text S4), is deterministic (no seed/MC noise),
  #   and resolves tail probabilities far below any Monte-Carlo floor.
  # method = "mc": the original Monte-Carlo Delta table (kept for comparison).
  #   NOTE the shipped MC also redrew the segment per state, computing an UNPAIRED
  #   statistic that inflates the FPR/FNR; that draw is now hoisted (see
  #   construct_Delta_table), so "mc" is the paired MC -- still noisy at small FPR.
  method <- match.arg(method)
  # Reproducibility: the only stochastic step is the Monte-Carlo Delta table
  # (rmultinom in construct_Delta_table). Seed it so a given object yields the
  # same rigidity every run; everything downstream (FPR/FNR, SE+/SE-, the
  # coarse->refine search) is deterministic given the Delta table. Restore the
  # caller's RNG state on exit so we don't perturb their stream. seed = NULL
  # leaves the RNG untouched (non-reproducible, old behaviour).
  if (!is.null(seed)) {
    old_seed = if (exists(".Random.seed", envir = .GlobalEnv))
      get(".Random.seed", envir = .GlobalEnv) else NULL
    set.seed(seed)
    on.exit(if (!is.null(old_seed)) assign(".Random.seed", old_seed, envir = .GlobalEnv), add = TRUE)
  }
  myDat = object
  seqlengths = seqlengths(myDat@Viterbi[[1]])
  expDesign = myDat@info$expDesign
  cat("file :", expDesign$files[1])
  if(is.null(average_coverage)){
    averages = c()
    for(samp in expDesign$files){
      f <- read.delim(file =samp, header = FALSE)
      averages = c(averages, mean(f$V4 + f$V6))
    }
    average_coverage = min(averages)
    # average_coverage = min(sapply(myDat@Viterbi, function(x) mean(x$total -1)))
  }
  if(is.null(crossovers_per_megabase)){
    crossovers_per_megabase = mean(sapply(myDat@Viterbi, function(x) length(rle(x$Viterbi)$values)/sum(seqlengths)))*1e6
    # crossovers_per_megabase = mean(colSums(calcCOnumber(myDat))/sum(seqlengths))
    crossovers_per_megabase =  max(calcCOnumber(myDat)/(seqlengths/1e6))
  }
  cat("Average coverage: ", average_coverage, "\n")
  cat("Crossovers/Mb: ", crossovers_per_megabase, "\n")
  # EXTRACT EMISSIONS
  CO_min = crossovers_per_megabase # the CO frequency used to construct the lower bound
  CO_max = 10*CO_min # the CO frequency used to construct the upper bound

  picked_parameters = myDat@params
  transition_pars = picked_parameters$transition  # this is on the absolute scale
  emission_pars = extract_emissions(picked_parameters)

  # Chromosome lengths for the genome length / chromosome count. Prefer the value
  # stored on @info (set by generateObject), but fall back to the lengths carried
  # by the Viterbi GRanges for objects built before that slot was populated.
  # Without this fallback, an empty slot gives total_length = 0, which zeroes every
  # segmentation-error term and makes optimize_R always return the smallest grid
  # rigidity regardless of the data.
  chromosome_lengths = myDat@info$seqlengths
  if (is.null(chromosome_lengths) || length(chromosome_lengths) == 0)
    chromosome_lengths = seqlengths
  n_chromosomes = length(chromosome_lengths)
  total_length = sum(chromosome_lengths)


# MAGIC CODE --------------------------------------------------------------
  states = c("mat","het","pat")

  # generate the grid for the u values on which FPR(u;r) and FNR(u;r) values will be evaluated
  max_segment_length = max_rigidity*2
  n_segments = 25
  segment_length_grid = pmin(2*ceiling(seq(1,(max_segment_length/2)^(1/4),length=n_segments)^4),max_segment_length)
  segment_length_grid = sort(unique(segment_length_grid))
  # ideally, we should use a grid with only ~3 significant binary digits (i.e., at most 3 columns in Delta_table
  # must be added up, this saves a lot of time)

  # Rigidity grid: powers of two from 2 up to max_rigidity. The manuscript sweeps
  # 2 <= r <= 400 (Text S4); the previous grid was 2^seq(from = 3, ...), i.e.
  # floored at 2^3 = 8, so it could never recommend r < 8 (and returned that floor
  # for low-coverage data whose true optimum is higher). A coarse power-of-two grid
  # for now; a later refinement can saturate finer points around the coarse minimum.
  rigidity_grid = as.integer(2^(1:floor(log2(max_rigidity))))
  n_rigidity = length(rigidity_grid)


  # generate the grid for the u values on which FPR(u;r) and FNR(u;r) values will be evaluated
  max_segment_length = max_rigidity*2
  n_segments = 25
  segment_length_grid = pmin(2*ceiling(seq(1,(max_segment_length/2)^(1/4),length=n_segments)^4),max_segment_length)
  segment_length_grid = sort(unique(segment_length_grid))
  # ideally, we should use a grid with only ~3 significant binary digits (i.e., at most 3 columns in Delta_table
  # must be added up, this saves a lot of time)


  # FPRfun(u, r) / FNRfun(u, r): the per-(segment length, rigidity) error rates,
  # supplied either by the exact tilted-FFT computation or the Monte-Carlo table.
  if (method == "exact") {
    em    <- emission_grid(average_coverage, emission_pars)
    logT  <- log(transition_pars)
    FPRfun <- function(u, r) exactFPR(u, r, em, logT)
    FNRfun <- function(u, r) exactFNR(u, r, em, logT)
  } else {
    Deltas = construct_Delta_table(n_obs = n_obs, # Monte Carlo draws for the FPR/FNR estimate (manuscript Text S4 eq. 53 uses N = 10^4 = the default; the original call wrongly passed the *biological sample* count, ~3, so the rates were estimated from 3 draws. Raise n_obs when the FPR is tiny -- e.g. high coverage / confident emissions -- and 1e4 draws leave SE_total noisy and the suggested r seed-unstable.)
                                   max_segment_length = max_segment_length, # length of the largest segment to be evaluated
                                   coverage=average_coverage, # average number of observations per marker
                                   emissions=emission_pars # named (mat,het,pat) list with emission probabilities (named vector c(alpha=...,beta=...) in each list entry)
    )
    FPRfun <- function(u, r) FPR(u, r, Deltas, transition_pars)$FPR
    FNRfun <- function(u, r) FNR(u, r, Deltas, transition_pars)$FNR
  }



  ### Segmentation-error model on a given rigidity grid.
  # eval_grid() computes the FPR/FNR grids and the SE+ / SE- / SE_total bounds for
  # an arbitrary set of rigidity values `rg`, reusing the single (expensive)
  # Monte-Carlo Delta table built above. Given `Deltas`, this is cheap and
  # deterministic, so we can evaluate as many rigidity values as we like -- which
  # is what makes the two-stage (coarse grid -> dense refine) search below free.
  eval_grid = function(rg){
    nr = length(rg)
    FPRg = matrix(0, n_segments, nr,
                  dimnames = list(paste0("Seg_", segment_length_grid), paste0("r_", rg)))
    FNRg = matrix(0, n_segments, nr,
                  dimnames = list(paste0("Seg_", segment_length_grid), paste0("r_", rg)))
    for (j in 1:nr) for (u in 1:n_segments){
      FPRg[u,j] = FPRfun(segment_length_grid[u], rg[j])
      FNRg[u,j] = FNRfun(segment_length_grid[u], rg[j])
    }
    # SE+ bounds (manuscript Text S4): the expected number of false segments is
    # bounded via a max over SEGMENT LENGTHS l (eq. 815: m = argmax_l FPR(l)/l;
    # eq. 830 lower: j = argmax_l FPR(l)). We take that max over segment_length_grid
    # and reuse FPRg -- NOT over the rigidity grid. The shipped code maxed over the
    # rigidity grid (a conflation of the l- and r-axes), which both made SE+ depend
    # on the grid's composition (so densifying it for the refine step distorted the
    # curve) and cost O(n_r^2) FPR evaluations. Using segment_length_grid makes SE+
    # a clean function of r, so the coarse->refine search is well-defined and cheap.
    up_SEp = sapply(1:nr, function(j) total_length * max(FPRg[,j] / segment_length_grid))
    lo_SEp = sapply(1:nr, function(j){
      jm = which.max(FPRg[,j])
      floor(total_length / segment_length_grid[jm]) * FPRg[jm,j]
    })
    lo_SEm = SEminus(total_length, n_chromosomes, COs_per_megabase = CO_min, FNRg, rg, segment_length_grid)
    up_SEm = SEminus(total_length, n_chromosomes, COs_per_megabase = CO_max, FNRg, rg, segment_length_grid)
    list(FPR_grid = FPRg, FNR_grid = FNRg,
         upper_SEplus = up_SEp, lower_SEplus = lo_SEp,
         upper_SEminus = as.numeric(up_SEm), lower_SEminus = as.numeric(lo_SEm),
         upper_SEtotal = up_SEp + as.numeric(up_SEm),
         lower_SEtotal = lo_SEp + as.numeric(lo_SEm))
  }

  ## Stage 1 -- coarse power-of-two grid.
  coarse_grid = rigidity_grid
  coarse = eval_grid(coarse_grid)
  ut = coarse$upper_SEtotal
  ic = which.min(ut); r_c = coarse_grid[ic]

  ## Refine only when the coarse optimum is interior AND the coarse SE_total curve
  ## is V-shaped (non-increasing up to the min, non-decreasing after). For a
  ## unimodal curve the continuous optimum is bracketed by the coarse minimum's two
  ## neighbours, [r_c/2, r_c*2].
  is_interior = ic > 1 && ic < length(coarse_grid)
  is_unimodal = all(diff(ut[1:ic]) <= 0) && all(diff(ut[ic:length(ut)]) >= 0)
  if (!is_interior){
    warning("optimize_R: coarse rigidity optimum is at the ",
            if (ic == 1) "lower" else "upper", " end of the grid (r = ", r_c,
            "); the true optimum may lie beyond it",
            if (ic == length(coarse_grid)) " -- consider increasing max_rigidity" else "",
            ". Returning the boundary value without refinement.", call. = FALSE)
  } else if (!is_unimodal){
    warning("optimize_R: coarse SE_total curve is not unimodal (Monte-Carlo noise?);",
            " returning the coarse optimum r = ", r_c, " without neighbour-bracket refinement.",
            call. = FALSE)
  }

  if (is_interior && is_unimodal){
    ## Stage 2 -- dense integer sweep over the neighbour bracket, reusing the same
    ## Delta table. Capped to <= 40 points so a wide bracket stays cheap. The coarse
    ## points are kept in the grid so SE+'s inner max still spans the full rigidity
    ## range. (A further iteration could re-bracket around this minimum; one pass
    ## suffices to pin an integer optimum here.)
    lo = max(2L, as.integer(r_c / 2))
    hi = min(as.integer(max_rigidity), as.integer(r_c * 2))
    refine = if (hi - lo + 1L <= 40L) lo:hi else unique(as.integer(round(seq(lo, hi, length.out = 40))))
    rigidity_grid = sort(unique(c(coarse_grid, as.integer(refine))))
    final = eval_grid(rigidity_grid)
  } else {
    rigidity_grid = coarse_grid
    final = coarse
  }

  ## Adopt the final grid's quantities (used by the suggestion, the return value,
  ## and the optional diagnostic plots below).
  n_rigidity   = length(rigidity_grid)
  FPR_grid     = final$FPR_grid
  FNR_grid     = final$FNR_grid
  upper_SEplus = final$upper_SEplus
  lower_SEplus = final$lower_SEplus
  lower_SEminus = final$lower_SEminus
  upper_SEminus = final$upper_SEminus
  upper_SEtotal = final$upper_SEtotal
  lower_SEtotal = final$lower_SEtotal
  indx = which.min(upper_SEtotal)
  rigidity_suggestion = rigidity_grid[indx]
  best_SEtotal = upper_SEtotal[indx]


# Plotting results (save_it = TRUE) ---------------------------------------

  if(save_it){
    options(scipen=999)
    n_rigidity = length(rigidity_grid)
    appendix = "Optimization-step"

    # FPR plot

    pdf(file.path(savedir, paste("FPR_plot_",appendix,".pdf",sep="")))
    plot(c(0,0), ylim=c(0,max(FPR_grid*100)),xlim=c(0.95,max(segment_length_grid)*2),
         type="n", xlab="Segment_length",ylab="FPR [%]",log="x")
    abline(h=0,col="grey")
    title(paste(c("FPR (segment_length,rigidity) , coverage = ",average_coverage),collapse=""))
    colorpal = rainbow(n_rigidity+3)[1:n_rigidity]
    lwidths = rep(2,n_rigidity) # (rigidity_grid==150)*1.5 + 1
    for (j in 1:n_rigidity){
      points(segment_length_grid,FPR_grid[,j]*100,type="l",col=colorpal[j],lwd = lwidths[j])
    }
    legend("topright",bty="n",legend = c("rigidity",rigidity_grid),col=c("white",colorpal),lty=1,lwd=2)
    dev.off()


    # FNR plot

    pdf(file.path(savedir, paste("FNR_plot_",appendix,".pdf",sep="")))
    plot(c(0,0), ylim=c(0,max(FNR_grid*100)),xlim=c(0.95,max(segment_length_grid)*2),
         type="n", xlab="Segment_length",ylab="FNR [%]",log="x")
    abline(h=0,col="grey")
    title(paste(c("FNR (segment_length,rigidity) , coverage = ",average_coverage),collapse=""))
    colorpal = rainbow(n_rigidity+3)[1:n_rigidity]
    lwidths = rep(2,n_rigidity) # (rigidity_grid==150)*1.5 + 1
    for (j in 1:n_rigidity){
      points(segment_length_grid,FNR_grid[,j]*100,type="l",col=colorpal[j],lwd = lwidths[j])
    }
    legend("topright",bty="n",legend = c("rigidity",rigidity_grid),col=c("white",colorpal),lty=1,lwd=2)
    dev.off()



    # Segmentation error per sample: SE+, SE-, SE_total plot

    mindisplay = 10^-5
    transform = function(x){log10(x+mindisplay)}
    backtransform = function(x){10^x-mindisplay}

    y_minmax = transform(c(0,10^5)) # max(c(upper_SEplus,lower_SEplus,upper_SEminus,lower_SEminus))))
    axpos1 = pretty(y_minmax)
    axpos = transform( backtransform(axpos1) + mindisplay ) [-1]
    axlabs = as.character(signif(backtransform(axpos),digits=3))

    x_minmax = range(c(rigidity_grid,rigidity_grid))

    pdf(file.path(savedir, paste("Segmentation_error_plot_",appendix,".pdf",sep="")),width=8,height=6)
    par(mar=c(5, 5.5, 4, 2.5) + 0.1)
    plot(rigidity_grid,transform(upper_SEplus), log="x",
         main=paste("Segmentation errors per sample\ncoverage = ",average_coverage,sep=""),
         xlab="Rigidity",
         ylab="",
         yaxt="n",
         ylim=y_minmax,
         xlim=x_minmax,
         type="n")
    title(ylab = paste("Wrong segments",sep=""),
          mgp = c(4, 3, 1))
    axis(side=2,at=axpos,labels=axlabs,las=1)
    axis(side=2,at=transform(0),labels=0,las=1)
    abline(h=transform(c(0,1)),col="grey")


    abline(v=rigidity_suggestion,lty=1,col="grey")

    linecolors = c(total="violet",pos="red",neg="blue")
    linewidths = c(total=3.5,pos=2,neg=1.8)

    # Total error (upper and lower bound)
    points(rigidity_grid,transform(upper_SEminus+upper_SEplus),type="l",
           col=linecolors["total"],lty=1,lwd=linewidths["total"])
    points(rigidity_grid,transform(lower_SEminus+lower_SEplus),type="l",
           col=linecolors["total"],lty=2,lwd=linewidths["total"])
    # False segment calls (upper and lower bound)
    points(c(rigidity_grid,x_minmax[2]),c(transform(upper_SEplus),transform(0)),type="l",
           col=linecolors["pos"],lwd=linewidths["pos"])
    points(rigidity_grid,transform(lower_SEplus),type="l",
           col=linecolors["pos"],lty=2,lwd=linewidths["pos"])
    # Missed segments (upper and lower bound)
    points(rigidity_grid,transform(upper_SEminus),type="l",
           col=linecolors["neg"],lty=1,lwd=linewidths["neg"])
    points(rigidity_grid,transform(lower_SEminus),type="l",
           col=linecolors["neg"],lty=2,lwd=linewidths["neg"])
    # add the text which marks the minimum error
    vpos = ifelse(best_SEtotal>1,10^-3,10)
    hpos = ifelse(rigidity_suggestion>10^4,2,4)
    text(rigidity_suggestion,transform(vpos),
         labels = paste("min total error at\nrigidity = ",rigidity_suggestion,sep=""),
         pos=4,offset=0.5)
    # add the legend for line colors and style
    legend("topright",legend=c("false segments","missed segements","total error",
                               "","upper bound","lower bound"),
           col = c(linecolors[c("pos","neg","total")],"white","dark grey","dark grey"),
           lty = c(1,1,1, 1, 1,2),lwd=2,bty="n",bg="white")
    dev.off()

    # save the results file
    resultsfile = file.path(savedir, paste("Simulation_results_",appendix,".RData",sep=""))
    save(
      average_coverage,
      rigidity_grid,
      upper_SEplus,
      lower_SEplus,
      upper_SEminus,
      lower_SEminus,
      upper_SEtotal,
      lower_SEtotal,
      segment_length_grid,
      FPR_grid,
      FNR_grid,
      rigidity_suggestion,
      best_SEtotal,
      file = resultsfile)
  }
  return(rigidity_suggestion)


}


#  Auxiliar functions -----------------------------------------------------

#'  utility function converting dec to reverse binary
#' @param num decimal numbers
#'
#' @keywords internal
#' @noRd
#'
#'

# utility function converting dec to reverse binary
dec2bin <- function(num){
  if (num %/% 2 == 0) return((num %% 2)==1)
  return(c((num %% 2==1),dec2bin(num %/% 2)))
} #dec2bin


#' Delta_table[x,y,m=segment_lengths,n_obs=#samples] contains, for each segment length and each observation, the probability for the observations in a given segment length generated by state x, evaluated as coming from y
#' @param n_obs HOw many samples shall be constructed
#' @param max_segment_length length of the largest segment to be evaluated
#' @param coverage average number of observations per marker
#' @param emissions named (mat,het,pat) list with emission probabilities (named vector c(alpha=...,beta=...) in each list entry)
#'
#' @keywords internal
#' @noRd
#'
# Delta_table[x,y,m=segment_lengths,n_obs=#samples] contains, for each segment length and each observation,
# the probability for the observations in a given segment length generated by state x, evaluated as coming from y
construct_Delta_table = function(n_obs = 10^4, # how many samples shall be constructed?
                                 max_segment_length = 2^10-1, # length of the largest segment to be evaluated
                                 coverage=1, # average number of observations per marker
                                 emissions # named (mat,het,pat) list with emission probabilities (named vector c(alpha=...,beta=...) in each list entry)
){

  # pre-calculate a lookup table <probs> of dimension (states,0:max_n_counts+1,0:max_n_counts+1)
  # with the entries probs[state,k,n] = Betabin(k-1;n-1,alpha[state],beta[state])

  states = c("mat","het","pat")
  max_n = ceiling(coverage * 2.5 + 8) # the largest number of observations per marker we cover (should never be seen)

  probs = array(NA,dim=c(3,max_n+1,max_n+1))
  dimnames(probs) = list(states,NULL,NULL)
  log_probs = array(NA,dim=c(3,max_n+1,max_n+1))
  dimnames(log_probs) = list(states,NULL,NULL)

  for (state in states){
    for (n in 0:max_n){
      probs[state,,n+1] = dbbinom(0:max_n, size=n,
                                  alpha= emissions[[state]]["alpha"],
                                  beta = emissions[[state]]["beta"], log=FALSE) * dpois(n,lambda=coverage)
      log_probs[state,,n+1] = dbbinom(0:max_n, size=n,
                                      alpha= emissions[[state]]["alpha"],
                                      beta = emissions[[state]]["beta"], log=TRUE) + dpois(n,lambda=coverage,log=T)
    } # end for n
  } # end for state
  log_probs[log_probs== -Inf] = 0

  # Initialize the Delta_table array
  max_m = floor(log2(max_segment_length))+1
  Delta_table = array(NA, dim = c(3,3,max_m,n_obs))
  dimnames(Delta_table) = list(states,states,NULL,NULL)

  # construct the Delta_table values for segments of length 2^(m-1) , m = 1,...,max_m
  #cat("Constructing Delta table with ",n_obs," observations.\n")
  #cat("Segment lengths (",max_m,"): ",sep="")
  for (m in 1:max_m){

    segment_length = 2^(m-1)
    #cat(segment_length,", ")

    for (x_state in states){
      # Draw the synthetic segments ONCE per (x_state, m) and score the SAME
      # observations under every y_state. The FPR/FNR are built from differences
      # like Delta[x,y] - Delta[x,x], which the manuscript defines as the PAIRED
      # log-likelihood ratio Delta_{m,x,y} = logP(o, all-y) - logP(o, all-x) on the
      # SAME observation o (Supplemental Text S4, eq. 60). Previously `tables` was
      # redrawn inside the y_state loop, so the two likelihoods were evaluated on
      # DIFFERENT random segments -> an UNPAIRED statistic whose variance is
      # inflated by +2*Cov(logP_y, logP_x). That overestimates FPR/FNR (e.g. an
      # FPR of ~1e-4 where the paired value is ~1e-10) and biases the suggested
      # rigidity upward. Drawing once and reusing it makes the estimate paired,
      # matching the manuscript's definition.
      tables = rmultinom(n_obs, size = segment_length, prob = probs[x_state,,])
      for (y_state in states){
        # evaluate these observations with the probabilities of the y_state
        Delta_table[x_state,y_state,m,] = colSums(tables * as.vector(log_probs[y_state,,]))
      } # end for y_state
    } # end for x_state

  } # end for m
  #cat("\n")

  return(Delta_table)
} # end construct_Delta_table

# ---- Exact paired FPR/FNR via exponential-tilted FFT convolution -------------
# These replace the Monte-Carlo Delta table when optimize_R(method = "exact").
# The FPR/FNR are tail probabilities P(D > t) of a sum D of per-marker
# log-likelihood-ratio increments (the manuscript's Delta_{m,x,y}, Supp Text S4),
# computed analytically: the increment lives on the tiny (k,n) support, so D's
# distribution is its L-fold self-convolution -- evaluated exactly by FFT. For
# the far tails that matter (FPR can be < 1e-15 at high coverage / confident
# emissions, unreachable by Monte-Carlo), the increment is exponentially TILTED
# by theta so the tilted L-fold mean equals the threshold; the convolution is then
# well-resolved near t and untilted via f_L(s) = g_L(s) * exp(L*K(theta) - theta*s).

#' Per-marker emission grid for the exact FPR/FNR (same probs/log_probs as
#' construct_Delta_table). `emissions` is the named (mat,het,pat) list.
#' @keywords internal
#' @noRd
emission_grid = function(coverage, emissions, states = c("mat","het","pat")){
  max_n = ceiling(coverage * 2.5 + 8)
  probs = array(NA, dim = c(3, max_n+1, max_n+1), dimnames = list(states, NULL, NULL))
  lp    = probs
  for (s in states) for (n in 0:max_n){
    probs[s,,n+1] = dbbinom(0:max_n, size = n, alpha = emissions[[s]]["alpha"],
                            beta = emissions[[s]]["beta"]) * dpois(n, lambda = coverage)
    lp[s,,n+1]    = dbbinom(0:max_n, size = n, alpha = emissions[[s]]["alpha"],
                            beta = emissions[[s]]["beta"], log = TRUE) + dpois(n, lambda = coverage, log = TRUE)
  }
  lp[lp == -Inf] = 0
  list(probs = probs, lp = lp)
}

#' Tail P(S {>=|>} thr) for S = sum over independent `groups`, each group a sum of
#' g$L iid per-marker increments g$x with weights g$w. Exponential tilting +
#' FFT convolution + untilt resolves tails far below any Monte-Carlo or plain-FFT
#' floor; deterministic. `strict` = TRUE gives P(S > thr), FALSE gives P(S >= thr).
#' @keywords internal
#' @noRd
tilt_tail = function(groups, thr, strict, H = 0.05){
  ats  = lapply(groups, function(g){ k = g$w > 0; list(w = g$w[k]/sum(g$w[k]), x = g$x[k], L = g$L) })
  Ltot = sum(vapply(ats, function(a) a$L, numeric(1)))
  xmin = Ltot * min(vapply(ats, function(a) min(a$x), numeric(1)))
  xmax = Ltot * max(vapply(ats, function(a) max(a$x), numeric(1)))
  if (thr > xmax) return(0)
  if (thr < xmin) return(1)
  lse    = function(v){ m = max(v); m + log(sum(exp(v - m))) }
  Ksum_p = function(th) sum(vapply(ats, function(a){ lw = log(a$w) + th*a$x; a$L * sum(exp(lw - lse(lw)) * a$x) }, numeric(1)))
  # saddlepoint: tilt so the tilted mean of S equals the threshold
  th = tryCatch(stats::uniroot(function(t) Ksum_p(t) - thr, c(-60, 60))$root, error = function(e) 0)
  allx = unlist(lapply(ats, `[[`, "x"))
  x0   = floor(min(allx)/H)*H
  M    = as.integer(round((max(allx) - x0)/H)) + 1L
  G    = 2L^ceiling(log2(Ltot*(M-1) + 2))
  cf   = rep(1+0i, G); logZ = 0
  for (a in ats){
    lw = log(a$w) + th*a$x; Kth = lse(lw); logZ = logZ + a$L*Kth
    h = numeric(G); idx = as.integer(round((a$x - x0)/H)) + 1L
    for (i in seq_along(idx)) h[idx[i]] = h[idx[i]] + exp(lw[i] - Kth)   # tilted, normalised
    cf = cf * fft(h)^a$L
  }
  gS   = Re(fft(cf, inverse = TRUE)) / G          # tilted L-fold convolution (centred at thr)
  vals = Ltot*x0 + (0:(G-1))*H
  # untilt and sum the tail in log space (only over resolved positive mass) so the
  # exp(logZ - theta*vals) factor never multiplies an FFT zero into a NaN.
  sel  = (if (strict) vals > thr else vals >= thr) & is.finite(gS) & gS > 0
  if (!any(sel)) return(0)
  s = sum(exp(log(gS[sel]) + logZ - th*vals[sel]))
  if (!is.finite(s)) return(0)
  min(1, max(0, s))
}

#' Exact paired FPR(u; r): probability a length-u flanking region is falsely
#' broken by inserting a central segment. em = emission_grid(); logT = log(transition).
#' @keywords internal
#' @noRd
exactFPR = function(u, r, em, logT, states = c("mat","het","pat")){
  if (u < r) return(0)
  m = matrix(0, 3, 3, dimnames = list(states, states))
  for (fl in states) for (ce in setdiff(states, fl)){
    const = logT[fl,ce] + logT[ce,fl] - 2*logT[fl,fl] + (u - r)*(logT[fl,fl] - logT[ce,ce])
    x = as.vector(em$lp[ce,,]) - as.vector(em$lp[fl,,])   # per-marker LLR increment (paired)
    m[fl,ce] = tilt_tail(list(list(w = as.vector(em$probs[fl,,]), x = x, L = u)), -const, TRUE)
  }
  0.25*sum(m["mat",]) + 0.25*sum(m["pat",]) + 0.5*sum(m["het",])
}

#' Exact paired FNR(u; r): probability a true central segment of length u is missed.
#' For u < r the statistic is a sum of two independent pieces (the length-u true
#' segment plus the length-(r-u) flanking remainder), handled as two groups.
#' @keywords internal
#' @noRd
exactFNR = function(u, r, em, logT, states = c("mat","het","pat")){
  m = matrix(NA_real_, 3, 3, dimnames = list(states, states))
  for (fl in states) for (ce in setdiff(states, fl)){
    const = logT[fl,ce] + logT[ce,fl] - 2*logT[fl,fl]
    xfc   = as.vector(em$lp[fl,,]) - as.vector(em$lp[ce,,])
    if (u >= r){
      thr = const + (u - r)*(logT[fl,fl] - logT[ce,ce])
      m[fl,ce] = tilt_tail(list(list(w = as.vector(em$probs[ce,,]), x = xfc, L = u)), thr, FALSE)
    } else {
      g = list(list(w = as.vector(em$probs[ce,,]), x = xfc, L = u),
               list(w = as.vector(em$probs[fl,,]), x = xfc, L = r - u))
      m[fl,ce] = tilt_tail(g, const, FALSE)
    }
  }
  (1/3)*m["mat","het"] + (1/3)*m["pat","het"] + (1/6)*sum(m["het", c("pat","mat")])
}

#' Compute the False positiv rate values
#' @param segment_length length of the segment for which the FPR is to be calculated
#' @param rigidity rigidity value of the rHMM
#' @param Delta_table the Delta table constructed by the construct_Delta_table funciton
#' @param transitions the transition matrix
#'
#' @keywords internal
#' @noRd
#'

FPR = function(segment_length, # length of the segment for which the FPR is to be calculated
               rigidity, # rigidity value of the rHMM
               Delta_table, # the Delta table constructed by the function above
               transitions # the transition matrix
){

  log_transitions = log(transitions)
  states = c("mat","het","pat")
  n_obs = dim(Delta_table)[4]
  FPR_xy = matrix(0,nrow=3,ncol=3,dimnames=list(states,states))
  if (segment_length<rigidity) return(list(FPR=0,FPR_xy=FPR_xy))

  # sum up the appropriate columns of the Delta_table values
  pick_m = which(dec2bin(segment_length))
  Delta_reduced = apply(Delta_table[,,pick_m,,drop=F],c(1,2,4),sum)

  for (flanking_state in states){
    for (central_state in setdiff(states,flanking_state)){

      # calculate the transition penalty for switching to another state
      constant = log_transitions[flanking_state,central_state] +
        log_transitions[central_state,flanking_state] -
        2 * log_transitions[flanking_state,flanking_state] +
        (segment_length-rigidity)*(log_transitions[flanking_state,flanking_state]-
                                     log_transitions[central_state,central_state])

      # calculate the relative frequency of the Delta value defined in the lyx being positive
      FPR_xy[flanking_state,central_state] = sum(Delta_reduced[flanking_state,central_state,] + constant >
                                                   Delta_reduced[flanking_state,flanking_state,]) / n_obs
    } # end central_state
  } # end flanking_state

  # weight the state-specific FPR rates to obtain a global FPR
  FPR = 1/4 * sum(FPR_xy["mat",]) + 1/4 * sum(FPR_xy["pat",]) + 1/2 * sum(FPR_xy["het",])
  return(list(FPR=FPR,FPR_xy=FPR_xy))
} # end FPR

#' Compte the False Negative Rate (FNR)
#' @param segment_length length of the segment for which the FPR is to be calculated
#' @param rigidity rigidity value of the rHMM
#' @param Delta_table the Delta table constructed by the construct_Delta_table funciton
#' @param transitions the transition matrix
#'
#' @keywords internal
#' @noRd
#'

FNR = function(segment_length, # length of the segment for which the FNR is to be calculated
               rigidity, # rigidity value of the rHMM
               Delta_table, # the Delta table constructed by the function above
               transitions # the transition matrix
){

  log_transitions = log(transitions)
  states = c("mat","het","pat")
  n_obs = dim(Delta_table)[4]
  FNR_xy = matrix(NA,nrow=3,ncol=3,dimnames=list(states,states))

  if (segment_length >= rigidity){

    pick_m = which(dec2bin(segment_length))
    Delta_reduced = apply(Delta_table[,,pick_m,,drop=F],c(1,2,4),sum)


    for (flanking_state in states){
      for (central_state in setdiff(states,flanking_state)){

        # calculate the transition penalty for switching to another state
        constant = log_transitions[flanking_state,central_state] +
          log_transitions[central_state,flanking_state] -
          2 * log_transitions[flanking_state,flanking_state] +
          (segment_length-rigidity)*(log_transitions[flanking_state,flanking_state]-
                                       log_transitions[central_state,central_state])

        # calculate the relative frequency of the Delta value defined in the lyx being positive
        FNR_xy[flanking_state,central_state] = sum(
          Delta_reduced[central_state,central_state,] + constant <=
            Delta_reduced[central_state,flanking_state,] ) / n_obs
      } # end central_state
    } # end flanking_state

  } else {
    # in case segment_length < rigidity do

    pick_m = which(dec2bin(segment_length))
    Deltay_reduced = apply(Delta_table[,,pick_m,,drop=F],c(1,2,4),sum)
    pick_mminusu = which(dec2bin(rigidity-segment_length))
    Deltax_reduced = apply(Delta_table[,,pick_mminusu,,drop=F],c(1,2,4),sum)

    for (flanking_state in states){
      for (central_state in setdiff(states,flanking_state)){

        # calculate the transition penalty for switching to another state
        constant = log_transitions[flanking_state,central_state] +
          log_transitions[central_state,flanking_state] -
          2 * log_transitions[flanking_state,flanking_state]

        # calculate the relative frequency of the Delta value defined in the lyx being positive
        FNR_xy[flanking_state,central_state] = sum(
          Deltay_reduced[central_state,central_state,] + Deltax_reduced[flanking_state,central_state,] + constant <=
            Deltay_reduced[central_state,flanking_state,] + Deltax_reduced[flanking_state,flanking_state,]) / n_obs
      } # end central_state
    } # end flanking_state
  } # end (if segment_length < rigidity)

  # weight the state-specific FNR rates to obtain a global FNR
  FNR = 1/3 * sum(FNR_xy["mat","het"]) + 1/3 * sum(FNR_xy["pat","het"]) + 1/6 * sum(FNR_xy["het",c("pat","mat")])
  return(list(FNR=FNR,FNR_xy=FNR_xy))
} # end FNR

#' Computations of the segment errors
#' @param total_length the total size of the genome (in bp)
#' @param n_chromosomes The number of chromosomes the organism has
#' @param COs_per_megabase The expected total number of COs per megabase
#' @param FNR_grid a matrix containing information of how many segments (rows) have been missed for each R (columns).
#' @param rigidity_grid The values of R for which it will be inspected the performance.
#'
#' @keywords internal
#' @noRd
#'

SEminus = function(total_length, # the total size of the genome (in bp)
                   n_chromosomes, # the number of chromosomes
                   COs_per_megabase, # the expected total number of COs per megabase
                   FNR_grid,
                   rigidity_grid,
                   segment_length_grid
){
  expected_COs = total_length/10^6*COs_per_megabase
  max_k = max(qpois(10^-7,expected_COs,lower.tail=FALSE),2)
  n_rigidity = length(rigidity_grid)
  n_segments = length(segment_length_grid)

  p_ulk = matrix(0, nrow = n_segments, ncol = max_k+1)
  for (k  in 1:max_k){
    p_ulk[,k+1] = diff(pbeta(c(0.5,segment_length_grid)/total_length,1,k+n_chromosomes))
  } # end for k

  k_vec = dpois(0:max_k, lambda = expected_COs)

  u_mat =  FNR_grid * diff(c(0,segment_length_grid))

  SE_minus = t(u_mat) %*% p_ulk %*% k_vec
  return(SE_minus)
} # end SEminus


