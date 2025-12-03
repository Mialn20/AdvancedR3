


load_description_packages <- function(desc_path = here::here("DESCRIPTION")) {
  if (!file.exists(desc_path)) {
    stop("DESCRIPTION file not found at: ", desc_path)
  }

  # Read DESCRIPTION
  desc <- read.dcf(desc_path)


  # Extract fields
  imports  <- if ("Imports"  %in% colnames(desc)) desc[1, "Imports"]  else ""
  suggests <- if ("Suggests" %in% colnames(desc)) desc[1, "Suggests"] else ""

  # Parse into package names
  pkg_list <- c(imports, suggests)
  pkg_list <- unlist(strsplit(pkg_list, "[,\n]+"))
  pkg_list <- trimws(pkg_list)
  pkg_list <- pkg_list[nzchar(pkg_list)]  # remove empty

  message("Found packages: ", paste(pkg_list, collapse = ", "))

  # Load all packages
  for (pkg in pkg_list) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      message(sprintf("Package '%s' is NOT installed.", pkg))
    } else {
      library(pkg, character.only = TRUE)
      message(sprintf("Loaded '%s'.", pkg))
    }
  }
}

lipidprocessesF <- function(x){
  return(x$gender)
}
