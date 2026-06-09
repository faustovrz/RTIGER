
#' Call Julia code to fit the values
#' @param rtigerobj an RTIGER object.
#' @param max.iter maximum number of iterations to acomplish by the EM.
#' @param eps differnece threshold to halt the EM.
#' @param trace logical value whether to trace the changes in the parameters along the iterations.
#' @param all logical value whether to use all data to fit the model.
#' @param random if all FALSE use random samples.
#' @param specific if all FALSE use specific samples.
#' @param nsamples if random TRUE, how many samples to use.
#' @param post.processing logical value, whether to run post.processing process.
#' @param progress_log optional file path. When non-NULL, Julia appends one
#'   newline-terminated record per EM iteration (iter/max, delta, eps, elapsed,
#'   per_iter, ETA<=) to this file, flushed each iteration so it can be tailed
#'   live for an ETA. NULL (default) disables it (Julia receives "" and writes
#'   nothing) — behavior is then identical to before. This is logging only and
#'   does not affect the fit. Note ETA<= is an upper bound (the EM usually
#'   converges at delta<eps before max.iter). Distinct from the developer
#'   debugInfo.txt dump.
#' @param verbose logical. When TRUE and progress_log is set, echo the progress
#'   records to the console after the fit returns (the fit call is blocking, so
#'   for a live view tail the file from a separate shell). Default FALSE.
#'
#' @return RTIGER object
#' @usage fit(rtigerobj, max.iter , eps,
#' trace, all = TRUE, random = FALSE,
#' specific = FALSE, nsamples = 20,
#' post.processing = TRUE, progress_log = NULL, verbose = FALSE)
#'
#' @examples
#'\dontrun{
#'data("fittedExample")
#' sourceJulia()
#' myfit = fit(myDat, max.iter = 2, eps=0.01,
#'             trace = TRUE, all = TRUE,
#'             random = FALSE, specific = FALSE,
#'             nsamples = 20, post.processing = TRUE)
#'
#'}
#' @export fit
#'

fit = function(rtigerobj, max.iter, eps, trace, all = TRUE, random = FALSE, specific = FALSE, nsamples = 20, post.processing = TRUE, progress_log = NULL, verbose = FALSE){
  params = rtigerobj@params
  obs = rtigerobj@matobs
  info = rtigerobj@info

  obs = lapply(obs, function(samp){
    chr = lapply(samp, function(cr){
      return(t(cr))
    })
    return(chr)
  })
  # Single on/off switch for the per-iteration progress log, orthogonal to
  # `verbose`: a non-empty path turns it on, "" keeps the fit silent (identical
  # to before). R owns the path. Accepts NULL/FALSE (off), a string (that path),
  # or TRUE (sentinel -> a default file under tempdir() when called directly;
  # RTIGER() resolves TRUE to <outputdir>/fit_progress.log before calling fit).
  progress_log_path =
    if (is.null(progress_log) || isFALSE(progress_log)) {
      ""
    } else if (isTRUE(progress_log)) {
      file.path(tempdir(), "fit_progress.log")
    } else {
      as.character(progress_log)
    }
  # cat("Inside fit the postprocessing value is:", post.processing, "\n")
  # function fit(Observations,info,initial_parameter,max_iter=100,eps=10^(-5),trace=false)
  myfit = julia_call("fit",obs, info, params, as.integer(max.iter), eps, trace, all, random , as.integer(nsamples), specific, post.processing, progress_log = progress_log_path )
  # function fit(
  #   input_Observations,
  #   info,
  #   initial_parameter,
  #   max_iter = 100,
  #   eps = 10^(-5),
  #   trace = false,
  #   all = true,
  #   random = true,
  #   nsamples=20,
  #   specific = nothing,
  # )
  nstates = myfit$parameterSet$nstates

  myvit = myfit$viterbiPath
  myvit = myvit[names(obs)]
  myvit = lapply(myvit, function(x) x[info$part_names])
  myvit = lapply(myvit, unlist)
  if(nstates == 3){
    vits = c( "pat", "het", "mat")
  } else if(nstates == 2){
    vits = c( 0,0)
    a = myfit$parameterSet$paraBetaAlpha
    b = myfit$parameterSet$paraBetaBeta
    par.dif = as.vector(a-b)
    hetst = which.min(abs(par.dif))
    vits[hetst] = "het"
    vits[-hetst] = ifelse(par.dif[-hetst] > 0, "pat", "mat")
  } else{
    vits = 1:nstates
  }


  for(i in info$sample_names){
    rtigerobj@Viterbi[[i]]$Viterbi = vits[myvit[[i]]]
  }
  rtigerobj@params = myfit$parameterSet
  rownames(rtigerobj@params$logtransition) = colnames(rtigerobj@params$logtransition) = vits
  rownames(rtigerobj@params$paraBetaAlpha) = vits
  rownames(rtigerobj@params$paraBetaBeta) = vits
  rownames(rtigerobj@params$logpi) = vits
  rownames(rtigerobj@params$pi) = vits
  rownames(rtigerobj@params$transition) = colnames(rtigerobj@params$transition) = vits

  rtigerobj@Probabilities = myfit[c("alpha", "beta", "gamma", "psi")]
  rtigerobj@num.iter = myfit$numberofiterations

  # verbose console echo of the progress records (the julia_call above blocks,
  # so this prints after completion; for a live ETA, tail the file meanwhile).
  if (verbose && nzchar(progress_log_path) && file.exists(progress_log_path)) {
    cat("EM progress (", progress_log_path, "):\n", sep = "")
    cat(readLines(progress_log_path), sep = "\n")
    cat("\n")
  }
  return(rtigerobj)

}
