#' frame_the_signal
#'
#' Realiza a divisão do sinal em short-frames
#'
#' @importFrom  purrr map
#' @return retorna um dataframe, onde cada linha representa o trecho do sinal dividido
#' @param signal um vetor numérico, representando o sinal a ser dividido
#' @param sample.rate representa o sample rate do sinal
#' @param frame.size representa o comprimento da janela em segundos.
#' @param frame.stride representa o passo de cada janela
#' @export
#' @examples
#' x <- sample(1000)
#' frames <- frame_the_signal(x, 1600)
#' head(frames)
#'
frame_the_signal <- function(signal, sample.rate, frame.size = 0.025,
                             frame.stride = 0.01) {
  # Transformar de segundos para amostras
  frame.length <- round(frame.size * sample.rate)
  frame.step <- round(frame.stride * sample.rate)
  signal.length <- length(signal)

  frames.amount <- ceiling( (signal.length - frame.length) / frame.step)

  # preencher com zeros o final do vetor para garantir que
  # os quadros possuam o mesmo comprimento
  pad.length <- frames.amount * frame.step + frame.length
  pad.signal <- c(signal, vector(mode = "numeric", pad.length -
                                   signal.length))

  # dividir o vetor em partes
  indices <- seq(from = 1, to = frames.amount * frame.step, by = frame.step)
  indices <- indices %>% purrr::map(function(x) x:(x + frame.length - 1))
  frames <- indices %>% purrr::map(function(x) pad.signal[x])

  do.call(rbind, frames)  #converter em um data.frame
  # frames
}


#' apply_window_hamming
#'
#' Realiza a aplicação da função Window Hamming que é dada pela seguinte equação:
#'
#' \deqn{w[n] = 0.56 - 0.46 * cos ((2\pi n)/N-1)}, onde o N representa o comprimento da janela.
#' @importFrom  signal hamming
#' @return retorna a lista de frames com a aplicação do window hamming
#' @param frames estrutura contendo os frames
#' @export
apply_window_hamming <- function(frames) {
  frame.length <- ncol(frames)
  windows <- apply(frames, 1, function(x) x * signal::hamming(frame.length))
  t(windows)
}


#' compute_power_spectrum
#'
#' Realiza a aplicação da função de Power Spectrum que é dada pela seguinte fórmula:
#'
#' \deqn{(|FFT(X_i)|^2)/N}, onde o N representa o comprimento da janela.
#'
#' @return retorna a lista de frames com a aplicação do Power Spectrum
#' @param frames estrutura contendo os frames
#' @param n valor representando o número de pontos para o cálculo da FFT
#' @export
compute_power_spectrum <- function(frames, n = 512) {
  P <- function(x) abs(nfft(x, n)) ^ 2 / n
  power.spectrum <- t(apply(frames, 1, P))
  power.spectrum
}

#' compute_mel_filterbanks
#'
#' Realiza o cálculo dos mel-filterbank.
#'
#' @importFrom  purrr map_dbl
#' @return retorna a lista de frames com a aplicação do Power Spectrum
#' @param freq.lower valor mínimo da frequência dada em hertz
#' @param freq.upper valor maximo da frequência dada em hertz
#' @param num.filters número de filtros considerados para o cálculo
#' @param npoints número de pontos considerado para o cálculo dos N-Points FFT
#' @param sample.rate sample rate do sinal
#' @export
compute_mel_filterbanks <- function(freq.lower, freq.upper, num.filters,
                                    npoints, sample.rate) {
  freq.lower <- convert_hertz_to_mel(freq.lower)
  freq.upper <- convert_hertz_to_mel(freq.upper)
  mel.points <- seq(from = freq.lower, to = freq.upper,
                    length.out = num.filters + 2)
  hertz.points <- mel.points %>% purrr::map_dbl(convert_mel_to_hertz)
  # To convert the frequencies to fft bin numbers we need to know
  # the FFT size and the sample rate
  fft.bin <- floor( (npoints + 1) * hertz.points / sample.rate)
  filter.banks <- 2:(num.filters + 1) %>%
    purrr::map(calcule_filter_bank,
               f = fft.bin, npoints = npoints)
  filter.banks <- do.call(rbind, filter.banks)
  filter.banks
}


#' mfcc_function
#'
#' Realiza o cálculo das MFFCs para um sinal.
#'
#' @return retorna a lista de frames com a aplicação do Power Spectrum
#' @param signal vetor numérico representando o sinal
#' @param sample.rate sample rate do sinal
#' @param freq.lower valor mínimo da frequência dada em hertz
#' @param freq.upper valor maximo da frequência dada em hertz
#' @param frame.size representa o comprimento da janela em segundos.
#' @param frame.stride representa o passo de cada janela
#' @param nfft número de pontos considerado para o cálculo dos N-Points FFT
#' @param num.filters número de filtros considerados para o cálculo
#' @param num.ceps representa o número de coeficientes a serem mantidos
#' @export
mfcc_function <- function(signal, sample.rate, freq.lower, freq.upper,
                          frame.size = 0.025, frame.stride = 0.01, nfft = 512,
                          num.filters = 40, num.ceps = 13) {
  # Aplicar filtro de pre-ênfase no sinal, valores padrões
  # de alfa são 0,95 ou 0,97
  emphasized_signal <- apply_preemphasis(signal, 0.97)

  # Passos do algoritmo:

  # 1.Dividir o sinal em short frames.
  frames <- frame_the_signal(emphasized_signal, sample.rate,
                             frame.size, frame.stride)
  #aplicar a função de hamming para cada frame
  frames <- apply_window_hamming(frames)

  # 2.Para cada quadro, calcular o power spectrum
  power.frames <- compute_power_spectrum(frames, n = nfft)

  # 3.Aplicar o mel filterbank aos power spectra,
  # somar a energia em cada filtro
  filter.banks <- compute_mel_filterbanks(freq.lower, freq.upper,
                                          num.filters, nfft, sample.rate)
  # Para calcular a energia do filter bank, multiplica-se
  # cada filter bank com seus power spectrum.
  filter.banks <- power.frames %*% t(filter.banks)

  # 4. Obter o logaritmo de todas as filterbank energies
  #substituir os zeros para evitar problemas com log
  filter.banks[filter.banks == 0] <- .Machine$double.eps
  filter.banks <- 20 * log10(filter.banks)

  # 5. Obter a DCT do log das filterbank energies.
  mfcc.var <- t(apply(filter.banks, 1, function(x) apply_dct(x)))

  # 6. Manter os coeficientes DCT 2-13, descartar o resto..
  mfcc.var <- mfcc.var[, 2:num.ceps]

  # Aplicar o sinusoidal liftering aos MFCCs
  apply_lifter(mfcc.var)
}
