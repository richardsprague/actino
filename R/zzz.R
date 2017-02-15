# zzz.R
# top level initializations

.onLoad <- function(libname = find.package("actino"), pkgname = "actino"){
  if(getRversion() >= "2.15.1")
    utils::globalVariables(c("ssr","tax_name","count","count_norm","tax_rank"))


  invisible()

}
