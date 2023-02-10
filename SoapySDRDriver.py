from pylab import psd
from pylab import show
import SoapySDR

from SoapySDR import * #SOAPY_SDR_ constants
import numpy #use numpy for buffers

#device = 'rtlsdr'
#Fs = 250e3
#frf = 915e6
def init_radio(device, Fs, frf, gain= None):
    
    args = dict(driver=device)
    sdr = SoapySDR.Device(args)

    #query device info
    print(sdr.listGains(SOAPY_SDR_RX, 0))
    if gain != None:
        print(sdr.setGain(SOAPY_SDR_RX, gain))

    #apply settings
    sdr.setSampleRate(SOAPY_SDR_RX, 0, Fs)
    sdr.setFrequency(SOAPY_SDR_RX, 0, frf)

    return sdr

#bufsize = 1024
def start_rf_stream(sdr, bufsize):

    #setup a stream (complex floats)
    rxStream = sdr.setupStream(SOAPY_SDR_RX, SOAPY_SDR_CF32)
    sdr.activateStream(rxStream) #start streaming

    #create a re-usable buffer for rx samples
    buff = numpy.array([0]*bufsize, numpy.complex64)

    return rxStream, buff

def get_radio_samples(sdr, rxStream, buff, printstats = True):
    sr = sdr.readStream(rxStream, [buff], len(buff))
    if printstats:
        print(sr.ret) #num samples or error code
        print(sr.flags) #flags set by receive operation
        print(sr.timeNs) #timestamp for receive buffer

    return buff

def stop_rf_stream(sdr, rxStream):
    sdr.deactivateStream(rxStream) #stop streaming
    sdr.closeStream(rxStream)

#Fs = 250000
#frf = 915e6
#NFFT = 1024
def simple_psd(buff, NFFT, Fs, frf):
    psd(buff, #data buffer
    NFFT, #FFT Size
    Fs/1e6, #Sampling Rate
    frf/1e6 #Center Frequency (Hz)
    ) 
    show()

'''
#example use
sdr = init_radio(device = 'rtlsdr', Fs = 250000, frf = 915e6, gain = None)
rxStream, buff = start_rf_stream(sdr, bufsize = 1024)

for i in range(5):
    get_radio_samples(sdr, rxStream, buff, printstats = False)

buff = get_radio_samples(sdr, rxStream, buff, printstats = False)

simple_psd(buff, NFFT = 1024, Fs = 250000, frf = 915e6)

stop_rf_stream(sdr, rxStream)
'''