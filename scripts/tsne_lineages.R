
check_envs <- function(x) {

  # make sure envs are shared between all files
  the_names <- sapply(x, function(x) sort(colnames(x)), simplify = FALSE)

  # test out
  intersect_mat <- sapply(the_names,
                          function(x) sapply(the_names,
                                             function(y) length(intersect(x,y))))
  # zero out diag elements
  intersect_mat <- intersect_mat * (1 - diag(dim(intersect_mat)[1]))

  # sum off-diagonal.
  off_diag_sum = sum(intersect_mat, na.rm = T)

  # check for correct sum
  target_sum <- max(sapply(the_names, length)) * length(the_names)

  if (off_diag_sum == target_sum) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}


do_cluster <- function(dat, meta) {

  set.seed(12345)

  tsne.norm = Rtsne(dat, pca = FALSE)

  return(tsne.norm)
}


grab_metadata <- function(x,
                          bcs,
                          target_names = c('Subpool.Environment',
                                           'Which.Subpools',
                                           'neutral_set')) {

  x2 <- x[x$Full.BC %in% bcs, names(x) %in% c('Full.BC', target_names)]

  # parse to source_env and ploidy
  x2$ploidy     <- sapply(x2$Subpool.Environment,
                          function(x) ifelse(grepl('2N',x),'2N','1N'))
  x2$source_env <- sapply(x2$Subpool.Environment,
                          function(s) {
                            tmp <- unlist(strsplit(s, '_'))
                            paste0(tmp[1:(length(tmp) - 1)], collapse = '_')
                          })

  x2$source_env <- sapply(x2$source_env,
                          function(x) ifelse(grepl('YPD', x), 'YPD', x))

  return(x2)
}
