#' Encode DNA Sequences as Categorical Nucleotide Positions
#'
#' Convert equal-length DNA sequences into a data frame of nucleotide
#' positions for use by the prediction model.
#'
#' @param dna_strings A non-empty character vector of equal-length DNA
#'   sequences containing uppercase A, T, C, and G, with no missing values.
#'   The supplied prediction model requires five-base sequences.
#'
#' @return A data frame with one row per sequence and columns named
#'   nt_pos1, nt_pos2, and so on. Each column is a factor with levels
#'   A, T, C, G, in that order.
#'
#' @keywords internal
#'
#' @examples
#' m6APrediction:::dna_encoding(c("GGACA", "AGACT"))
dna_encoding <- function(dna_strings) {
  nn <- nchar(dna_strings[1])

  seq_m <- matrix(
    unlist(strsplit(dna_strings, "")),
    ncol = nn,
    byrow = TRUE
  )

  colnames(seq_m) <- paste0("nt_pos", 1:nn)
  seq_df <- as.data.frame(seq_m)

  seq_df[] <- lapply(
    seq_df,
    factor,
    levels = c("A", "T", "C", "G")
  )

  return(seq_df)
}


#' Predict m6A Probabilities and Statuses for Multiple Sites
#'
#' Apply a fitted random forest classifier to a table of site features.
#' DNA sequences are encoded as categorical nucleotide positions, and
#' RNA categories are assigned the factor levels used during training.
#'
#' @param ml_fit A fitted randomForest classification model with a
#'   class named Positive. The package includes a compatible model.
#' @param feature_df A non-empty data frame containing the following
#'   columns, with no missing values:
#'   \describe{
#'     \item{gc_content}{Numeric GC content between 0 and 1.}
#'     \item{RNA_type}{One of mRNA, lincRNA, lncRNA, or pseudogene.}
#'     \item{RNA_region}{One of CDS, intron, 3'UTR, or 5'UTR.}
#'     \item{exon_length}{Numeric exon-length feature on the log2 scale
#'       used in the course dataset.}
#'     \item{distance_to_junction}{Numeric junction-distance feature on
#'       the log2 scale used in the course dataset.}
#'     \item{evolutionary_conservation}{Numeric conservation score
#'       between 0 and 1.}
#'     \item{DNA_5mer}{Five uppercase DNA bases using A, T, C, and G.
#'       The supplied model was trained on contexts with central A.}
#'   }
#' @param positive_threshold A single numeric probability cutoff between
#'   0 and 1. The default is 0.5. Scores strictly greater than this
#'   cutoff are Positive; scores equal to or below it are Negative.
#'
#' @return The input data frame with RNA_type and RNA_region converted
#'   to factors, plus two columns: predicted_m6A_prob, a numeric
#'   positive-class probability, and predicted_m6A_status, a character
#'   value of Positive or Negative.
#'
#' @details Probabilities are obtained from the random forest's
#'   positive-class vote proportions. Predictions reflect learned
#'   associations and do not establish experimental confirmation.
#'   The bundled example table demonstrates software use and is not
#'   an independent dataset for evaluating model accuracy.
#'
#' @import randomForest
#' @importFrom stats predict
#' @export
#'
#' @examples
#' ml_fit <- readRDS(system.file(
#'   "extdata", "rf_fit.rds",
#'   package = "m6APrediction", mustWork = TRUE
#' ))
#' example_df <- read.csv(system.file(
#'   "extdata", "m6A_input_example.csv",
#'   package = "m6APrediction", mustWork = TRUE
#' ))
#' prediction_multiple(
#'   ml_fit,
#'   example_df[1:3, ],
#'   positive_threshold = 0.6
#' )
prediction_multiple <- function(
    ml_fit,
    feature_df,
    positive_threshold = 0.5
) {
  stopifnot(all(
    c(
      "gc_content",
      "RNA_type",
      "RNA_region",
      "exon_length",
      "distance_to_junction",
      "evolutionary_conservation",
      "DNA_5mer"
    ) %in% colnames(feature_df)
  ))

  feature_df$RNA_type <- factor(
    feature_df$RNA_type,
    levels = c("mRNA", "lincRNA", "lncRNA", "pseudogene")
  )

  feature_df$RNA_region <- factor(
    feature_df$RNA_region,
    levels = c("CDS", "intron", "3'UTR", "5'UTR")
  )

  inferred_prob <- predict(
    ml_fit,
    newdata = cbind(
      feature_df,
      dna_encoding(feature_df$DNA_5mer)
    ),
    type = "prob"
  )

  feature_df$predicted_m6A_prob <- inferred_prob[, "Positive"]
  feature_df$predicted_m6A_status <- "Negative"

  feature_df$predicted_m6A_status[
    feature_df$predicted_m6A_prob > positive_threshold
  ] <- "Positive"

  return(feature_df)
}


#' Predict m6A Probability and Status for a Single Site
#'
#' Construct a one-row feature table and call prediction_multiple()
#' to obtain the positive-class probability and predicted status.
#'
#' @param ml_fit A fitted randomForest classification model with a
#'   class named Positive.
#' @param gc_content A single numeric GC content between 0 and 1.
#' @param RNA_type A single RNA category: mRNA, lincRNA, lncRNA,
#'   or pseudogene.
#' @param RNA_region A single RNA region: CDS, intron, 3'UTR,
#'   or 5'UTR.
#' @param exon_length A single numeric exon-length feature on the
#'   log2 scale used in the course dataset.
#' @param distance_to_junction A single numeric junction-distance
#'   feature on the log2 scale used in the course dataset.
#' @param evolutionary_conservation A single numeric conservation
#'   score between 0 and 1.
#' @param DNA_5mer A single five-base uppercase DNA sequence using
#'   A, T, C, and G. The supplied model was trained on contexts
#'   with A in the central position.
#' @param positive_threshold A single numeric cutoff between 0 and 1.
#'   Scores strictly greater than the cutoff are Positive.
#'   Equality is classified as Negative. The default is 0.5.
#'
#' @return A named character vector of length two with elements
#'   predicted_m6A_prob and predicted_m6A_status. The probability is
#'   converted to text because both values share one vector.
#'   Use as.numeric() on the probability element for calculations.
#'
#' @seealso \code{\link{prediction_multiple}}
#' @export
#'
#' @examples
#' ml_fit <- readRDS(system.file(
#'   "extdata", "rf_fit.rds",
#'   package = "m6APrediction", mustWork = TRUE
#' ))
#' prediction_single(
#'   ml_fit,
#'   gc_content = 0.6,
#'   RNA_type = "mRNA",
#'   RNA_region = "CDS",
#'   exon_length = 12,
#'   distance_to_junction = 5,
#'   evolutionary_conservation = 0.8,
#'   DNA_5mer = "GGACA",
#'   positive_threshold = 0.5
#' )
prediction_single <- function(
    ml_fit,
    gc_content,
    RNA_type,
    RNA_region,
    exon_length,
    distance_to_junction,
    evolutionary_conservation,
    DNA_5mer,
    positive_threshold = 0.5
) {
  feature_df <- data.frame(
    gc_content = gc_content,
    RNA_type = factor(
      RNA_type,
      levels = c("mRNA", "lincRNA", "lncRNA", "pseudogene")
    ),
    RNA_region = factor(
      RNA_region,
      levels = c("CDS", "intron", "3'UTR", "5'UTR")
    ),
    exon_length = exon_length,
    distance_to_junction = distance_to_junction,
    evolutionary_conservation = evolutionary_conservation,
    DNA_5mer = DNA_5mer
  )

  pred_outcome <- prediction_multiple(
    ml_fit = ml_fit,
    feature_df = feature_df,
    positive_threshold = positive_threshold
  )

  returned_vector <- c(
    pred_outcome$predicted_m6A_prob,
    pred_outcome$predicted_m6A_status
  )

  names(returned_vector) <- c(
    "predicted_m6A_prob",
    "predicted_m6A_status"
  )

  return(returned_vector)
}
