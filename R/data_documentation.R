#' Intestinal Parasites in European Passerine Birds
#'
#' A dataset examining the prevalence of intestinal parasites in the fecal droppings of
#' 366 European migratory and non-migratory passerine birds, originally published in Bandelj et al.
#'
#' The dataset includes information on species-level variation and phylogenetic classification,
#' along with covariates describing migratory behavior and dietary type. It has been used as
#' a motivating example for hierarchical logistic regression with random effects.
#'
#' @format A data frame with 366 rows and 6 variables:
#' \describe{
#'   \item{ID.bird}{Unique identifier for each bird (integer).}
#'   \item{parasites}{Binary outcome indicating presence (1) or absence (0) of intestinal parasites.}
#'   \item{Phylogenetic.Tomi}{Phylogenetic classification (factor) of the bird. There are 8 clusters ranging from 1 to 158 birds per cluster.}
#'   \item{migration}{Migratory behavior (factor) with levels "migratory" and "non-migratory". 298 migratory and 68 non-migratory birds.}
#'   \item{food}{Diet type (factor) with levels "granivorous", "mixed", "omnivorous", and "insectivorous". 13, 211, 25, and 117 birds, respectively.}
#'   \item{species}{Species identifier (factor), nested within the phylogenetic classification. 42 species clusters ranging from 1 to 56 birds per cluster.}
#' }
#'
#' @details
#' The dataset can be used to fit generalized linear mixed models (GLMMs) for binary outcomes,
#' including random intercepts for species and both random intercepts and slopes for migration
#' at the phylogenetic level. This hierarchical structure allows for modeling of
#' phylogenetic and species-level dependencies while evaluating the effects of migration
#' and diet on parasite prevalence.
#'
#' @source
#' Bandelj et al., Original publication dataset on intestinal parasites in European passerine birds.
"birds"
