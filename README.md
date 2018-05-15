# Implementação dos MFCC em R

O presente repositório contém a implementação dos Mel-frequency cepstral coefficients (MFCCs) em R. O pacote `MFCC` foi desenvolvido utilizando o controle de versão do Github. Então você pode fazer download dele utilizando o seguinte comando no terminal do R:

``` r
# install.packages("devtools")
devtools::install_github("JessicaSousa/MFCC")
```

### Exemplo

Inicialmente, para carregar um arquivo de áudio é necessário ter instalado no R o pacote `tuneR`, esse pacote pode ser facilmente instalado através do comando `install.packages('tuneR')`. Com as bibliotecas corretamente configuradas, pode-se calcular os MFCCs usando a função `mfcc_function`, a qual recebe os seguintes parâmetros:

| Parâmetros        | Descrição           | 
| ------------- |:-------------:|   
|signal| O sinal de áudio o qual será realizada a extração das MFCCs|
|sample.rate| O sample rate do sinal de áudio|
|freq.lower| Menor valor de frequência dada em hertz|
|freq.upper| Maior valor frequência dada em hertz|
|frame.size| Representa o comprimento da janela, em segundos|
|frame.stride| Representa o passo de cada janela, em segundos|
|nfft| Número de pontos considerado para o cálculo dos N-Points FFT|
|num.filters| Quantidade de filtros considerados no filterbank|
|num.ceps| Representa até onde serão mantidos os coeficientes (2-13, se num.ceps = 13)|

``` r
#Carregar arquivo de áudio
#Carregar os 5 primeiros segundos do áudio
sound.data <- tuneR::readWave('audio.wav', from = 0, to = 5, units = "seconds") 
sound <- sound.data@left # sinal de áudio
sample.rate <- sound.data@samp.rate 

#Calcular as MFCCs para o sinal
mfcc <- MFCC::mfcc_function(sound.data, sample.rate, freq.lower = 0, 
                            freq.upper = sample.rate / 2,
                            frame.size = 0.025, frame.stride = 0.01,
                            nfft = 512, num.filters = 40, num.ceps = 13) 
```


#### Exemplos mais detalhados
Exemplos com gráficos presente no seguinte [link](https://jessicasousa.github.io/MFCC/inst/doc/example.html)

### Referências

[Practical Cryptography](http://practicalcryptography.com/miscellaneous/machine-learning/guide-mel-frequency-cepstral-coefficients-mfccs/)
[SR Wiki](http://recognize-speech.com/feature-extraction/mfcc)
