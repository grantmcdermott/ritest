.onAttach <- function(libname, pkgname){
  options("pboptions" = list(
    type = if (interactive()) "timer" else "none",
    char = "-",
    txt.width = 50,
    gui.width = 300,
    style = 3,
    initial = 0,
    title = "R progress bar",
    label = "",
    nout = 5L,
    min_time = 2))
  invisible(NULL)
}
