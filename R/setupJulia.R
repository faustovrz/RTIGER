#' Installs the needed packages in JULIA to run the EM algorithm for rHMM.
#'
#' @param JULIA_HOME the file folder which contains julia binary, if not set, JuliaCall will look at the global option JULIA_HOME, if the global option is not set, JuliaCall will then look at the environmental variable JULIA_HOME, if still not found, JuliaCall will try to use the julia in path.
#'
#' @return empty
#' @usage setupJulia(JULIA_HOME = NULL)
#'
#' @export setupJulia
setupJulia = function(JULIA_HOME = NULL){
  if(!is.null(JULIA_HOME)) julia_setup(JULIA_HOME = JULIA_HOME)
  v = julia_eval("string(VERSION)")
  # The optimize-julia-core fork is developed and validated on Julia 1.12.6
  # (native arm64 on Apple Silicon). Only warn for genuinely older Julia; newer
  # is fine (and faster). Proper component-wise compare via numeric_version, so
  # e.g. 1.9 is correctly treated as older than 1.12 (a naive 1.9 > 1.12 numeric
  # test would be wrong).
  rec = "1.12.6"
  vshort = sub("^([0-9]+\\.[0-9]+\\.[0-9]+).*$", "\\1", v)
  older = tryCatch(numeric_version(vshort) < numeric_version(rec),
                   error = function(e) FALSE)
  if (older)
    cat(sprintf("Note: Julia %s detected. RTIGER's optimized core is validated on Julia %s or newer (native arm64 on Apple Silicon); older Julia may be slower or hit incompatibilities.\n", v, rec))
  julia_install_package_if_needed("Optim")
  julia_install_package_if_needed("Distributions")
  julia_install_package_if_needed("LinearAlgebra")
  julia_install_package_if_needed("CSV")
  julia_install_package_if_needed("DelimitedFiles")
  julia_install_package_if_needed("DataFrames")
  # julia_install_package_if_needed("Plots")
}

#' Function needed before using RTIGER() function. It loads the scripts in Julia that fit the rHMM.
#'
#' @return empty
#' @export sourceJulia
#'

sourceJulia = function(){
  julia_source(paste(system.file("julia", package = "RTIGER"),"/AuxilaryFunctions.jl", sep = ""))
  julia_source(paste(system.file("julia", package = "RTIGER"), "/rHMM_methods.jl", sep =""))
  # julia_source(paste(system.file("julia", package = "RTIGER"), "/rHMM_methods_old.jl", sep =""))

}
