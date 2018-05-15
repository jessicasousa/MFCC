library(magrittr)

#' kronecker_delta
#'
#' A função kronecker delta e definida como:
#'
#' @return retorna um número de acordo com o valor i
#' @param i corresponde um valor inteiro.
#' @param j corresponde um valor inteiro.
#' @export
#' @examples
#' print(kronecker_delta(i = 1))
#'
kronecker_delta <- function(i, j = 0) {
    if (i != j)
        return(0)
    return(1)
}


#' dct1
#'
#' O cálculo da DCT-1, é implementado conforme consta na documentação do matlab.
#'
#' https://www.mathworks.com/help/signal/ref/dct.html
#' @importFrom  purrr map_dbl
#' @return retorna o cálculo da DCT-1 para um vetor numérico.
#' @param k representa o valor usado para o kronecker delta
#' @param x representa o vetor numérico
#' @param N representa o valor que corresponde ao comprimento do vetor x.
#' @export
dct1 <- function(k, x, N) {
    f <- function(n) {
        x[n] / sqrt(1 + kronecker_delta(n, 1) + kronecker_delta(n, N)) *
        cos( (pi * (n - 1) * (k - 1)) / (N - 1))
    }
    n <- 1:N
    cst <- sqrt(2 / N) * (1 / sqrt(1 + kronecker_delta(k, 1) +
                                 kronecker_delta(k, N)))
    cst * sum(n %>% purrr::map_dbl(f))
}

#' dct2
#'
#' O cálculo da DCT-2, é implementado conforme consta na documentação do matlab.
#'
#' https://www.mathworks.com/help/signal/ref/dct.html
#' @importFrom  purrr map_dbl
#' @return retorna o cálculo da DCT-2 para um vetor numérico.
#' @param k representa o valor usado para o kronecker delta
#' @param x representa o vetor numérico
#' @param N representa o valor que corresponde ao comprimento do vetor x.
#' @export
dct2 <- function(k, x, N) {
    f <- function(n) {
        x[n] * cos( (pi * (k - 1) * (2 * n - 1)) / (2 * N))
    }
    n <- 1:N
    cst <- sqrt(2 / N) * 1 / sqrt(1 + kronecker_delta(k, 1))
    cst * sum(n %>% purrr::map_dbl(f))
}

#' dct3
#'
#' O cálculo da DCT-3, é implementado conforme consta na documentação do matlab.
#'
#' https://www.mathworks.com/help/signal/ref/dct.html
#' @importFrom  purrr map_dbl
#' @return retorna o cálculo da DCT-3 para um vetor numérico.
#' @param k representa o valor usado para o kronecker delta
#' @param x representa o vetor numérico
#' @param N representa o valor que corresponde ao comprimento do vetor x.
#' @export
dct3 <- function(k, x, N) {
    f <- function(n) {
        x[n] / sqrt(1 + kronecker_delta(n, 1)) * cos( (pi * (2 * k - 1) *
                                                    (n - 1)) / (2 * N))
    }
    n <- 1:N
    cst <- sqrt(2 / N)
    cst * sum(n %>% purrr::map_dbl(f))
}

#' dct4
#'
#' O cálculo da DCT-4, é implementado conforme consta na documentação do matlab.
#'
#' https://www.mathworks.com/help/signal/ref/dct.html
#' @importFrom  purrr map_dbl
#' @return retorna o cálculo da DCT-4 para um vetor numérico.
#' @param k representa o valor usado para o kronecker delta
#' @param x representa o vetor numérico
#' @param N representa o valor que corresponde ao comprimento do vetor x.
#' @export
dct4 <- function(k, x, N) {
    f <- function(n) {
        x[n] * cos( (pi * (2 * k - 1) * (2 * n - 1)) / (4 * N))
    }
    n <- 1:N
    cst <- sqrt(2 / N)
    cst * sum(n %>% purrr::map_dbl(f))
}


#' apply_dct
#'
#'
#' Dado um vetor, realiza o cálculo da DCT, a implementação e baseada
#' nessa documentação do Matlab. No qual são utilizadas as DCT-I, DCT-II,
#' DCT-III e DCT-IV.
#'
#' https://www.mathworks.com/help/signal/ref/dct.html
#'
#' @importFrom  purrr map_dbl
#' @return retorna o sinal transformado pela DCT
#' @param x representa o sinal
#' @param type representa o tipo de transformação DCT utilizada (type = 'I',
#' type = 'II', type = 'III' ou type = 'IV').
#' @export
apply_dct <- function(x, type = "II") {
    N <- length(x)
    k <- 1:N
    switch(type,
           I = k %>% purrr::map_dbl(dct1, x = x, N = N),
           II = k %>% purrr::map_dbl(dct2, x = x, N = N),
           III = k %>% purrr::map_dbl(dct3, x = x, N = N),
           IV = k %>% purrr::map_dbl(dct4, x = x, N = N))
}
