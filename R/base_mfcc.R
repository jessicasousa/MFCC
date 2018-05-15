#' Pipe operator
#'
#' @name %>%
#' @rdname pipe
#' @keywords internal
#' @export
#' @importFrom magrittr %>%
#' @usage lhs \%>\% rhs
NULL


#' MFCC: Mel-Frequency Cepstral Coefficients.
#'
#' Trabalho da disciplina Fundamentos de Computação Gráfica, implementação de
#' Mel-Frequency Cepstral Coefficients (MFCCs).
#'
#' @docType package
#' @name MFCC
#' @author Jessica Cardoso
#' @references
#' \itemize{
#' \item https://www.mathworks.com/help/signal/ref/dct.html
#' \item http://practicalcryptography.com/miscellaneous/machine-learning/guide-mel-frequency-cepstral-coefficients-mfccs/
#' \item https://docs.scipy.org/doc/numpy-1.14.0/reference/routines.fft.html
#' }
NULL



#' convert_hertz_to_mel
#'
#'A escala Mel tem como objetivo imitar a percepção não-linear do ouvido humano, sendo
#'mais discriminativa nas freqüências mais baixas e menos discriminativa nas freqüências
#'mais altas. A conversão da frequência hertz para a escala mel é dada pela seguinte fórmula:
#' \deqn{Mel(hz) = 2595 * log10(1 + hz / 700)}
#'
#' @return Tem como retorno o valor em escala mel.
#' @param hz É um valor numérico que representa o valor da frequência em hertz
#' @export
#' @examples
#'
#' mel <- convert_hertz_to_mel(1600)
#' print(mel)
#'
convert_hertz_to_mel <- function(hz) {
    2595 * log10(1 + hz / 700)
}

#' convert_mel_to_hertz
#'
#' Realiza a conversão da escala mel para a frequência em hertz. A fórmula é dada por:
#' \deqn{Hz(mel) = 700 * (10 ^ (mel / 2595) - 1)}
#'
#' @return A função retorna o valor numérico da frequência em hertz
#' @param mel É um valor numérico que representa o valor da escala Mel
#' @export
#' @examples
#'
#' hz <- convert_hertz_to_mel(1600)
#' print(hz)
convert_mel_to_hertz <- function(mel) {
    700 * (10 ^ (mel / 2595) - 1)
}


#' apply_preemphasis
#'
#' O filtro de pré-ênfase é utilizado com o objetivo de amplificar o sinal de
#' altas frequências e é util de várias maneiras como:
#' \enumerate{
#' \item equilibrar o espectro de frequências, uma vez que as frequências altas costumam ter magnitudes
#' menores em comparação com frequências mais baixas,
#' \item evitar problemas numéricos durante a operação de transformada de Fourier
#' \item também pode melhorar a relação sinal-ruido (SNR - Signal-to-Noise Ratio ).
#' }
#' O filtro pode ser aplicado usando a seguinte equação:
#'
#' \deqn{y(t) = x(t) - \alpha x(t-1)}
#'
#'
#' os valores padrões do coeficiente \eqn{\alpha} são 0,95 ou 0,97
#'
#' @importFrom  purrr map_dbl
#' @return Retorna um vetor numérico representando o sinal após a aplicação do filtro
#' de pre-ênfase.
#' @param x refere-se a um vetor numérico representando o sinal de som.
#' @param a refere-se ao coeficiente do filtro \eqn{\alpha} , default, \eqn{\alpha}  = 0,95
#' @export
#' @examples
#' x <- sample(50) #signal
#' print(x)
#' apply_preemphasis(x)
apply_preemphasis <- function(x, a = 0.95) {
    indices <- 1:length(x)
    indices %>% purrr::map_dbl(function(t) if (t == 1)
      x[t] else x[t] - a * x[t - 1])
}


#' nfft
#'
#' Esta função realiza o cálculo da Transformada Discreta de Fourier (DFT) sobre
#' n-pontos de um vetor unidimensional utilizando o algoritmo da Transformada
#' Rapida de Fourier (FFT). A DFT sobre valores reais tem como resultado um sinal
#' Hermitian symmetry, ou seja, os termos da frequência negativa sao apenas conjugados
#' complexos dos termos de frequência positiva. Os termos de valores negativos sao
#' descartados do resultado do cálculo para essa função.
#'
#' Sobre DFT ser Hermitian symmetry:
#' https://docs.scipy.org/doc/numpy/reference/routines.fft.html
#' @importFrom  stats fft
#' @return Retorna um vetor de comprimento teto de n/2, representando a aplicação da
#' transformada de Fourier sobre o sinal de entrada.
#' @param signal refere-se a um vetor numérico representando o sinal o qual se deseja
#' aplicar a transformada.
#' @param n parâmetro opcional que representa a quantidade de pontos a ser considerado para o
#' cálculo da DFT. Se nenhum valor for repassado, o cálculo será realizado sobre todos os elementos
#' do sinal. Se \code{n} for maior que o comprimento, o sinal e preenchido com zeros, se for menor,
#' o sinal e cortado.
#' @export
#' @examples
#' x <- sample(50) #signal
#' print(x)
#' nfft(x)
nfft <- function(signal, n = Inf) {
    signal.length <- length(signal)
    if(n == Inf){
      n = signal.length
    }
    if (n > signal.length) {
        # pad signal with zeros
        signal <- c(signal, vector(mode = "numeric", n - signal.length))
    }
    if(signal.length )
    nfft_var <- fft(signal[1:n])
    n <- ceiling( (n + 1) / 2)
    nfft_var[1:n]
}


#' filter_bank_function
#'
#' Esta função realiza o cálculo sobre um valor numérico k de acordo com a fórmula definida na
#' seção \strong{Computing the Mel filterbank} do
#' seguinte \href{http://http://www.practicalcryptography.com/miscellaneous/machine-learning/guide-mel-frequency-cepstral-coefficients-mfccs//}{tutorial}.
#'
#' @return Retorna um valor numérico de acordo com o intervalo que se encontra o valor k.
#' @param k refere-se a um valor real.
#' @param f um vetor de números reais representando as frequências convertidas de números de
#' bin fft.
#' @param m consiste em um número inteiro que varia no intervalo de 1 até os números de pontos
#' calculados na FFT .
#' @export
filter_bank_function <- function(k, f, m) {
    if (k < f[m - 1] || k > f[m + 1]) {
        return(0)
    } else if (k >= f[m - 1] && k < f[m]) {
        return( (k - f[m - 1]) / (f[m] - f[m - 1]))
    }
    return( (f[m + 1] - k) / (f[m + 1] - f[m]))
}

#' calcule_filter_bank
#'
#' Esta função calcula o filter bank para um sinal aplicando um respectivo filtro \code{m}. Onde \code{m}
#' varia de 1 ate a quantidade filtros desejadas.
#'
#' @importFrom  purrr map_dbl
#' @return Retorna um vetor numérico representando o resultado do sinal sobre aquele filtro.
#' @param m refere-se a um número real representando o filtro, este valor varia de 1 ate a quantidade
#' de filtros desejadas.
#' @param f um vetor de números reais representando as frequências convertidas de números de
#' bin fft.
#' @param npoints consiste na quantidade de pontos calculados na FFT.
#' @export
calcule_filter_bank <- function(m, f, npoints) {
    indices <- (f[m - 1]:f[m + 1])
    Hm <- vector(mode = "numeric", npoints / 2 + 1)
    Hm[indices + 1] <- indices %>% purrr::map_dbl(filter_bank_function,
                                                  f = f, m = m)
    Hm
}

#' apply_lifter
#'
#' Esta função aplica o sinusoidal liftering aos MFCCs para realçar os MFCCs mais elevados.
#' E é definida como a seguinte equação:
#' \deqn{\hat{MFCC_i} = 1 + ((w_i * D)/2)sin((\pi * n)/D)}
#'
#' @return Retorna a mfcc com a aplicação do sinusoidal liftering.
#' @param mfcc uma matriz representando a MFCC, onde cada coluna representa os coeficientes da DCT.
#' @param cepstral.lifter consiste na quantidade de pontos calculados na FFT.
#' @export
apply_lifter <- function(mfcc, cepstral.lifter = 22) {
    n <- 0:(ncol(mfcc) - 1)
    wi <- 1 + (cepstral.lifter / 2) * sin(pi * n / cepstral.lifter)
    sweep(mfcc, MARGIN = 2, wi, `*`)
}
