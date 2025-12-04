


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


create_table_descriptive_stats <- function(x){
  x1 <- x |>
    dplyr::group_by(metabolite) |>
    dplyr::summarise(across(value, list(mean = mean, sd = sd))) |>
    dplyr::mutate(across(where(is.numeric), \(x) round(x, digits = 1))) |>
    dplyr::mutate(MeanSD = glue::glue("{value_mean} ({value_sd})")) |>
    dplyr::select(Metabolite = metabolite, `Mean SD` = MeanSD)
  return(x1)
}

create_plot_distributions <- function(data) {
  data |>
  ggplot2::ggplot(ggplot2::aes(x = value)) +
    ggplot2::geom_histogram() +
    ggplot2::facet_wrap(ggplot2::vars(metabolite), scales = "free") +
    ggplot2::theme_minimal()
}



ListProcess <- function(x) {
  # Example processing: calculate the sum of each vector
  result <- sum(x)
  return(result)
}


#' Do some cleaning to fix issues in the data.
#'
#' @param data
#'
#' @returns clean data frame
clean <- function(data) {
  data |>
    dplyr::group_by(dplyr::pick(-value)) |>
    dplyr::summarise(value = mean(value), .groups = "keep") |>
    dplyr::ungroup()
}

#' Fix data to process it for model fitting.
#'
#' @param data The lipidomics data.
#'
#' @returns A data frame.
#'
preprocess <- function(data) {
  data |>
    dplyr::mutate(
      class = as.factor(class),
      value = scale(value)
    )
}

#' fit model
#'
#' @param data
#' @param model_formula
#'
#' @returns
#' @export
#'
#' @examples
fit_model <- function(data,model_formula){
  glm(
    formula = model_formula,
    data = data,
    family = binomial
  ) |>
    broom::tidy(exponentiate = TRUE) |>
    dplyr::mutate(
      metabolite = unique(data$metabolite),
      model = format(model_formula),
      .before = tidyselect::everything()
    )
}

#' Title
#'
#' @param data
#'
#' @returns
#' @export
#'
#' @examples
create_model_results <- function(data){
  data |>
    dplyr::group_split(metabolite) |>
    purrr::map(preprocess) |>
    purrr::map(fit_all_models) |>
    purrr::list_rbind()
}


#' fit_all_models function to a given dataframe
#'
#' @param data
#'
#' @returns
#' @export
#'
#' @examples
fit_all_models <- function(data) {
  list(
    class ~ value,
    class ~ value + gender + age
  ) |>
    purrr::map(\(model) fit_model(data, model = model)) |>
    purrr::list_rbind()
}



#' Making a plot function
#'
#' @param results
#'
#' @returns
#' @export
#'
#' @examples
create_plot_model_results <- function(results){
  results |>
    dplyr::filter(term == "value",std.error <= 2, estimate < 5) |>
    dplyr::select(metabolite, model, estimate, std.error,p.value) |>
    ggplot2::ggplot(ggplot2::aes(x = estimate, y = metabolite,xmin=estimate - std.error, xmax = estimate + std.error))+
    ggplot2::geom_pointrange()+
    ggplot2::geom_vline(xintercept = 1, linetype = "dashed",color="red")+
    ggplot2::facet_grid(cols = ggplot2::vars(model))+
    ggplot2::theme_minimal()+
    ggplot2::labs(x="odds ratio", y="Metabolite")+
    ggplot2::geom_text(ggplot2::aes(label = paste("p.value: ",round(p.value,5),"(don't look here, Luke)"), nudge_y = 0.2))
}





