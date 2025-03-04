groups <- f$ls(recursive = TRUE)
for (name in groups$name[grepl("categories", groups$name)]) {
  names <- strsplit(name, "/")[[1]]
  names <- c(names[1:length(names) - 1], "levels")
  new_name <- paste(names, collapse = "/")
  f[[new_name]] <- f[[name]]
}
for (name in groups$name[grepl("codes", groups$name)]) {
  names <- strsplit(name, "/")[[1]]
  names <- c(names[1:length(names) - 1], "values")
  new_name <- paste(names, collapse = "/")
  f[[new_name]] <- f[[name]]
  grp <- f[[new_name]]
  grp$write(args = list(1:grp$dims), value = grp$read() + 1)
}
f$close_all()
