#!/usr/bin/python
#from scipy.io import wavfile
import scipy.io.wavfile as wavfile
#import numpy.fft.fftpack as fft
from numpy.fft import fft
from matplotlib.pyplot import *
from numpy import *

def savePlot(name, data):
	plot(data)
	savefig(name)
	cla()

def subset(x, i, N, padding='zeros'):
	assert padding in ('zeros', 'circular')
	assert N <= len(x)
	if padding == 'zeros':
		if i < N//2:
			x_subset = x[:i + N//2]
			x_subset = concatenate([zeros(N//2 - i), x_subset])
		else:
			x_subset = x[i - N//2 : i + N//2]
			if len(x_subset) < N:
				x_subset = concatenate([x_subset, zeros(N - len(x_subset))])
	elif padding == 'circular':
		temp = concatenate([x, x, x])
		center = len(x) + i
		x_subset = temp[center-N//2:center+N//2]
	return x_subset


def zeropad(w, N):
	T = len(w)
	w_zeros = zeros(N//2 - T//2)
	w = concatenate([w_zeros, w, w_zeros])
	if len(w) == N + 1:
		w = w[:-1]
	return w


def stft(x, w, L=None):
	# L is the overlap, see http://cnx.org/content/m10570/latest/
	N = len(x)
	#T = len(w)
	if L is None:
		L = N
	# Zerro pad the window
	w = zeropad(w, N)
	X_stft = []
	points = range(0, N, N//L)
#	points = arange(,1227+1)
	print len(points)
	
	for i in points:
		x_subset = subset(x, i, N)
		fft_subset = fft.fft(x_subset * w)
		X_stft.append(fft_subset)
	print "done"
	X_stft = array(X_stft).transpose()
	return X_stft


def spec(x, w):
	return abs(stft(x, w))


#1 - leer el archivo con las muestras
#	el resultado de read es una tupla, el elemento 1 tiene las muestras
the_file = wavfile.read("az.wav")
#the_file = wavfile.read("busytone.wav")
rate = the_file[0]
length = len(the_file[1])
data = the_file[1]

power = 10
while pow(2,power) < length:
	power += 1

data = append(data, zeros(pow(2,power) - length))

print len(data)

savePlot("original_signal.png",data)

#2 - aplico transformada de fourier
trans = fft.rfft(data)
trans_real = abs(trans)

#2b - lo grafico
savePlot("transformed.png",trans_real)

#3 - busco la frecuencia
band = 1000
frec = argmax(trans_real)
filter_array = append(zeros(frec - (band / 2)), ones(band))
filter_array = append(filter_array, zeros(len(trans_real) - len(filter_array)))

print frec
filtered_array = multiply(trans, filter_array)

savePlot("filtered_trans.png",abs(filtered_array))


#4 - antitransformo
filtered_signal = array(fft.irfft(filtered_array)[:length], dtype="int16")
savePlot("filtered_signal.png",filtered_signal)


#5 - grabo el archivo filtrado
wavfile.write("out.wav", rate,filtered_signal)

#6 ??
#spectrogram = spec(filtered_signal, ones(64))
#imshow(spectrogram)
#axis(xmin=0,xmax=2000)
spectrogram = specgram(filtered_signal)
savefig("spectogram.png")
cla()

print len(spectrogram[0])

# spectrogram[0] es la matriz del rojo
red_matrix = spectrogram[0]
vec_ones = ones(len(red_matrix))
vec_sum = (matrix(vec_ones) * matrix(red_matrix)).transpose()
savePlot("frecuency_volume.png",vec_sum)


presence = zeros(len(vec_sum))

threshold = max(vec_sum) / 2.0
for i in range(len(presence)):
	if vec_sum[i] > threshold:
		presence[i] = 1

plot(presence)
axis(ymin=0,ymax=10)
savefig("presence.png", dpi=300)
cla()

print len(vec_sum)
