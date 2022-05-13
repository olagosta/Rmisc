#' Binary similarity coefficients
#'
#' Computes similarity coefficient for two binary variables.
#'
#' @param i,j Numeric, character or logical vectors of same length with same
#'   two tokens. If a vector contains only one of the tokens, them must be
#'   type factor with both tokens as levels.
#' @param method Character string indicating which coefficient to be computed.
#'    Default to "phi". Can be abbreviated.
#' @param missing.value Value to be treated as NA.
#'
#' @details
#' Coefficients are calculated based on a 2x2 contingency table or confusion matrix:
#'
#'          j                     predicted
#'        1   0            a        P   N
#'       -------           c      ---------
#'   1  | a | b |    OR    t   P | TP | FN |
#' i     -------           u      ---------
#'   0  | c | d |          a   N | FP | TN |
#'       -------           l      ---------
#'
#' Where:
#' a = number of cases on which i and j are 1 (positive matches / true positive)
#' b = number of cases where i is 1 and j is 0 (mismatches - i occurrence / false negative)
#' c = number of cases where i is 0 and j is 1 (mismatches - j occurrence / false positive)
#' d = number of cases where both i and j are 0 (negative matches / true negative)
#' a + b + c + d = n, the number of variables
#'
#' Compreensive lists of similarity metrics can be found in:
#' https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0247751
#' http://www.iiisci.org/journal/PDV/sci/pdfs/GS315JG.pdf
#'
#' In cases where the contingency table have zeros, some coefficients will
#' return Inf or NaN. That's not an error, you just need a coefficient more
#' adequated for your data.
#'
#' @return A numeric value between 0 and 1 or -1 and 1, depending on the coefficient.
#'
#' @export
#'
#' @example
#' i <- c(0,0,0,1,1,1)
#' j <- c(0,1,1,1,0,0)
#' similarity(i, j, method = "maxwell")

similarity <- function(i, j,
  method = c("phi", "ari", "baroni-urbani", "consonni-todeschini", "dennis",
             "dispersion", "hamann", "maxwell-pilliner", "michael", "peirce1",
             "peirce2", "peirce3", "scott", "tarantula", "yuleQ", "yuleW"),
  missing.value = NA) {


    if (length(i) != length(j)) stop("i and j must have the same length")

    if (!is.na(missing.value)) {
      i[i == missing.value] <- NA
      j[j == missing.value] <- NA
    }

    if (!is.numeric(method)) method <- match.arg(method)


  tbl <- table(i,j)

    if (length(tbl) != 4) stop("i and j must have the same binary levels")

    a <- tbl[2,2]
    b <- tbl[2,1]
    c <- tbl[1,2]
    d <- tbl[1,1]
    n <- sum(tbl) # a+b+c+d


  # ARI
  if (method == "ari") return ((n*(a+d)-((a+b)*(a+c)+(c+d)*(b+d)))/(n^2-((a+b)*(a+c)+(c+d)*(b+d))))

  # Baroni-Urbani and Buser II
  if (method == "baroni-urbani-buser2") return ((sqrt(a*d)+a-b-c)/(sqrt(a*d)+a+b+c))

  # Consonni and Todeschini
  if (method == "consonni-todeschini") return ((log(1+a*d)-log(1+b*c))/log(1+n^2/4))

  # Dennis
  if (method == "dennis") return ((a*d-b*c)/sqrt(n*(a+b)*(a+c)))

  # Dispersion
  if (method == "dispersion") return ((a*d-b*c)/n^2)

  # Hamann
  if (method == "hamann") return ((a+d-b-c)/n)

  # Maxwell and Pilliner
  if (method == "maxwell-pilliner") return (2*(a*d-b*c)/((a+b)*(c+d)+(a+c)*(b+d)))

  # Michael
  if (method == "michael") return (4*(a*d-b*c)/((a+d)^2+(b+c)^2))

  # Peirce I
  if (method == "peirce1") return ((a*d-b*c)/((a+b)*(c+d)))

  # Peirce II
  if (method == "peirce2") return ((a*d-b*c)/((a+c)*(b+d)))

  # Peirce III (0:1)
  if (method == "peirce3") return ((a*b+b*c)/(a*b+2*b*c+c*d))

  # Phi / Matthews Correlation Coefficient (MCC)
  if (method == "phi" | method == "mcc") return ((a*d-b*c)/sqrt((a+b)*(a+c)*(b+d)*(c+d)))

  # Scott
  if (method == "scott") return ((4*a*d-(b+c)^2)/((2*a+b+c)*(2*d+b+c)))

  # Tarantula (0:1)
  if (method == "tarantula") return ((a*(c+d))/(c*(a+b)))

  # Yule’s Q
  if (method == "yuleQ") return ((a*d-b*c)/(a*d+b*c))

  # Yule’s W
  if (method == "yuleW") return ((sqrt(a*d)-sqrt(b*c))/(sqrt(a*d)+sqrt(b*c)))
}
