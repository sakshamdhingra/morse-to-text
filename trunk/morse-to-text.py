#!/usr/bin/python
import scipy.io.wavfile as wavfile
from numpy.fft import fft
from matplotlib.pyplot import *
from numpy import *

def savePlot(name, data):
	plot(data)
	savefig(name)
	cla()


class SoundFile:

	def __init__(self, path):
		#1 - leer el archivo con las muestras
		#	el resultado de read es una tupla, el elemento 1 tiene las muestras
		the_file = wavfile.read(path)
		self.rate = the_file[0]
		self.length = len(the_file[1])
		self.data = the_file[1]
		# appendea ceros hasta completar una potencia de 2
		power = 10
		while pow(2,power) < self.length:
			power += 1
		self.data = append(self.data, zeros(pow(2,power) - self.length))
		print len(self.data)
	
	def setdata(self, data):
		self.data = data

	def getdata(self):
		return self.data

	def getlength(self):
		return self.length

	def saveas(self, path):
		wavfile.write(path, self.rate, self.data)

	def savePlot(self, fileName):
		savePlot(fileName,data)


class SignalFilter:

	def filter(self, soundfile):
		#2 - aplico transformada de fourier
		trans = fft.rfft(soundfile.getdata())
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
		filtered_signal = array(fft.irfft(filtered_array)[:soundfile.getlength()], dtype="int16")
		savePlot("filtered_signal.png",filtered_signal)
		soundfile.setdata(filtered_signal)

class SpectreAnalyzer:

	def spectrogram(self, signal):
		spectrogram = specgram(signal)
		savefig("spectrogram.png")
		cla()
		print len(spectrogram[0])
		return spectrogram

	def sumarizecolumns(self, mat):
		vec_ones = ones(len(mat))
		vec_sum = (matrix(vec_ones) * matrix(mat)).transpose()
		savePlot("frecuency_volume.png",vec_sum)
		return vec_sum

	def findpresence(self, vec_sum):
		presence = zeros(len(vec_sum))
		threshold = max(vec_sum) / 2.0
		for i in range(len(presence)):
			if vec_sum[i] > threshold:
				presence[i] = 1
		plot(presence)
		axis(ymin=0,ymax=10)
		savefig("presence.png", dpi=300)
		cla()
		return presence

	def findpulses(self, soundfile):
		spec = self.spectrogram(soundfile.getdata())
		# spec[0] es la matriz del rojo
		red_matrix = spec[0]
		vec_sum = self.sumarizecolumns(red_matrix)
		presence = self.findpresence(vec_sum)
		return presence

the_file = SoundFile("az.wav")
the_filter = SignalFilter()
the_filter.filter(the_file)
the_file.saveas("filtered.wav")

analyzer = SpectreAnalyzer()
analyzer.findpulses(the_file)

