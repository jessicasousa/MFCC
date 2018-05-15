## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ------------------------------------------------------------------------
#install.packages('tuneR')
library(MFCC)
library(tuneR)

#Carregar arquivo de áudio exemplo
sound.data <- MFCC::sound.data
sound <- sound.data@left #valores do arquivo de áudio
sample.rate <- sound.data@samp.rate #sample rate do arquivo de áudio
sound.time <- 0:(length(sound)-1)/sample.rate

## ----fig.height = 2, fig.width = 7, fig.align = "center"-----------------
library(ggplot2)

data.raw <- data.frame(x = sound.time, y = sound)

ggplot(data.raw, aes(x, y)) + geom_line() +
  xlab("Tempo (s)") + ylab("Amplitude") #+


## ----fig.height = 2, fig.width = 7, fig.align = "center"-----------------

emphasized_signal <- apply_preemphasis(sound, 0.97)

data.emphasized <- data.frame(x = sound.time, y = emphasized_signal)

#Exibição
ggplot(data.emphasized, aes(x, y)) + geom_line() +
  xlab("Tempo (s)") + ylab("Amplitude")

## ------------------------------------------------------------------------
fft.npoints <- 512 #números de pontos considerados para o cálculo da fft 
#help(nfft) #para mais informações
freq.lower <- 0 #frequência mínima em hertz considerada
freq.upper <- sample.rate / 2 #frequência máxima em hertz considerada
num.filters <- 40 #número filtros considerados para o filterbank

#1.Dividir o sinal em short frames.
frames <- frame_the_signal(emphasized_signal, sample.rate)
frames <- apply_window_hamming(frames) #aplicar a função de hamming para cada frame

#2.Para cada quadro, calcular o power spectrum
power.frames <- compute_power_spectrum(frames, n = fft.npoints)

#3.Aplicar o mel filterbank aos power spectra, somar a energia em cada filtro
fbanks <- compute_mel_filterbanks(freq.lower,freq.upper, num.filters, fft.npoints, sample.rate)
#Para calcular a energia do filter bank, multiplica-se cada filter bank com seus power spectrum.
filter.banks <- power.frames %*% t(fbanks)

#4. Obter o logaritmo de todas as filterbank energies
filter.banks[filter.banks == 0] <- .Machine$double.eps #substituir os zeros para evitar problemas com log
filter.banks <- 20 * log10(filter.banks)

## ----fig.height = 2, fig.width = 7, fig.align = "center"-----------------
x <- seq(from = freq.lower, to = freq.upper, length.out = ncol(fbanks)) %>% rep(num.filters)
y <- t(fbanks)
y %<>% as.data.frame() %>% tidyr::gather()

data <- data.frame(x = x, values = y$value, filters = y$key)

ggplot(data, aes(x, values, colour = filters)) +
  geom_line() + xlab("Frequência") + ylab("Amplitude") + theme(legend.position="none") +
  scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0))

## ----fig.height = 2, fig.width = 7, fig.align = "center"-----------------

#Organizar dado para melhor visualização
fbanks.spec <- reshape2::melt(filter.banks)
fbanks.spec$Var1 <- fbanks.spec$Var1 / 100 
fbanks.spec$Var2 <- fbanks.spec$Var2 / 10 

ggplot(fbanks.spec, aes(Var1,Var2)) + geom_raster(aes(fill = value)) +
  scale_fill_gradientn(colours = rainbow(10)) +
  xlab("Tempo (s)") + ylab("Frequência (kHz)") + ggtitle("Espectograma do sinal") +
  scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0)) 


## ------------------------------------------------------------------------
#5. Obter a DCT do log das filterbank energies.
mfcc <- t(apply(filter.banks, 1, function(x) apply_dct(x)))

#6. Manter os coeficientes DCT 2-13, descartar o resto..
mfcc <- mfcc[, 2:13]

## ----fig.height = 2, fig.width = 7, fig.align = "center"-----------------
#Aplicar o sinusoidal liftering aos MFCCs
mfcc.lift <- apply_lifter(mfcc)


#Organizar dado para melhor visualização
mfccs.spec <- reshape2::melt(mfcc.lift)
mfccs.spec$Var1 <- mfccs.spec$Var1 / 100 

#Espectograma do MFCCs
ggplot(mfccs.spec, aes(Var1,Var2, fill=value)) + geom_raster(aes(fill = value)) +
  scale_fill_gradientn(colours = rainbow(10)) +
  xlab("Tempo (s)") + ylab("Coeficientes das MFCCs") + ggtitle("MFCCs") +
  scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0)) 



