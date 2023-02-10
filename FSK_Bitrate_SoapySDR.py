import numpy as np
import scipy.fft
import scipy.signal
import matplotlib.pyplot as plt
#from rtlsdr import RtlSdr
import SoapySDRDriver
from pylab import psd
from pylab import show


def moving_average(a, n) :
    ret = np.cumsum(a, dtype=float)
    ret[n:] = ret[n:] - ret[:-n]
    return ret[n - 1:] / n

def get_time_from_samples(arr, samplerate):
    if type(arr) != list:
        arrlength = arr.size
    else:
        arrlength = len(arr)

    return arrlength / samplerate

def specto(sig, bins, graph, fft_sizecust = None):
    testsignal = sig
    fft_size = 2**bins
    num_rows = int(np.floor(len(testsignal)/fft_size))
    if fft_sizecust != None:
        spectrogram = np.zeros((num_rows, fft_sizecust))
    else: 
        spectrogram = np.zeros((num_rows, fft_size))
    for i in range(num_rows):
        if fft_sizecust != None:
            spectrogram[i,:] = 10*np.log10(np.abs(np.fft.fftshift(np.fft.fft(testsignal[i*fft_size:(i+1)*fft_size], fft_sizecust)))**2)
        else:
            spectrogram[i,:] = 10*np.log10(np.abs(np.fft.fftshift(np.fft.fft(testsignal[i*fft_size:(i+1)*fft_size])))**2)
    if graph == True:
        if fft_sizecust != None:
            plt.imshow(spectrogram, aspect='auto', extent = [0, fft_sizecust, len(testsignal), 0])
        else:
            plt.imshow(spectrogram, aspect='auto', extent = [0, fft_size, len(testsignal), 0])
        plt.xlabel("Frequency [MHz]")
        plt.ylabel("Time [s]")
        plt.show()

    return spectrogram

def IQ_conj(IQarr):
    return np.conjugate(IQarr)

def IQ_delay(IQarr, delay, real = None):
    arrreal = IQarr.real
    arrimag = IQarr.imag
    #rolling by delay indexes
    arrrealdelay= np.roll(arrreal, delay)
    arrimagdelay = np.roll(arrimag, delay)

    if real:
        return arrrealdelay
    else:
        return arrrealdelay + 1j * arrimagdelay
    
#this function uses a quadrature demodulation algorithm to create an instantaneous frequency signal
def frequency_detector(samples, samplerate, averaging = None, graph = None):
    #delaying samples by 1
    samples_delay = IQ_delay(samples, 1)

    #getting complex conjugate signal
    samples_delay_conj = IQ_conj(samples_delay)

    #np multiply implements complex multiplication
    samples_mult = np.multiply(samples, samples_delay_conj)

    #taking the angle of the previous signal gives us instantaneous frequency (then multiplying by a constant to get the correct frequency) 
    samples_angle = np.angle(samples_mult)
    freq = np.multiply((samplerate / 6.3), samples_angle)

    if averaging != None:
        freq = moving_average(freq, averaging)

    if graph != None:
        plt.plot(freq)
        plt.show()

    return freq

#makes an IQ signal using hilbert transform. currently not used
def make_IQ_signal(sig):
    sigimag = scipy.signal.hilbert(sig)
    return sig + 1j * sigimag

#detects signal based on fft power from spectrogram. Currently not used.
def detectsignal(sm, airtime, Fs, resolution, fftpow):
    spectrogram = specto(sm, resolution, graph = None)

    fft_size = 2**resolution
    num_rows = int(np.floor(len(sm)/fft_size))
    
    ilst = []
    binlst = []
    for i in range(num_rows):
        rowmax = np.max(spectrogram[i,:])
        print(rowmax)
        if rowmax > fftpow:
            ilst.append(i)
            symbin = np.argmax(spectrogram[i,:])
            binlst.append(symbin)

    try:
        smsig = sm[ilst[0]*fft_size:ilst[-1]*fft_size]
        sigairtime = get_time_from_samples(smsig, Fs)
        if sigairtime > airtime:
            return smsig  
    except:
        pass

    return np.array([])

#gets average based on max and min
def get_avg(sig, avgnum):
    siglen = sig.size
    xinit = 3000
    sigmax = np.array([])
    sigmin = np.array([])
    x = xinit
    while x < siglen - avgnum - xinit:
        sigarr = sig[x:x+avgnum]
        sigmax = np.append(sigmax, np.max(sigarr))
        sigmin = np.append(sigmin, np.min(sigarr))

        x += avgnum
    return (np.mean(sigmax) + np.mean(sigmin)) / 2

#next three functions are methods to capture clock phase/frequency for coherent recieving.    
def midpoint(a):
    mean_a = np.mean(a)
    mean_a_greater = np.ma.masked_greater(a, mean_a)
    high = np.ma.median(mean_a_greater)
    mean_a_less_or_equal = np.ma.masked_array(a, ~mean_a_greater.mask)
    low = np.ma.median(mean_a_less_or_equal)
    return (high + low) / 2

def find_clock_frequency(spectrum):
    maxima = scipy.signal.argrelextrema(spectrum, np.greater_equal)[0]
    while maxima[0] < 2:
        maxima = maxima[1:]
    if maxima.any():
        threshold = max(spectrum[2:-1])*0.8
        indices_above_threshold = np.argwhere(spectrum[maxima] > threshold)
        return maxima[indices_above_threshold[0]]
    else:
        return 0

def wpcr(a):
    if len(a) < 4:
        return []
    b = (a > midpoint(a)) * 1.0
    d = np.diff(b)**2
    if len(np.argwhere(d > 0)) < 2:
        return []
    f = scipy.fft.fft(d, len(a))
    p = find_clock_frequency(abs(f))
    if p == 0:
        return []
    cycles_per_sample = (p*1.0)/len(f)
    clock_phase = 0.5 + np.angle(f[p])/(np.pi * 2)
    if clock_phase <= 0.5:
        clock_phase += 1
    symbols = []
    for i in range(len(a)):
        if clock_phase >= 1:
            clock_phase -= 1
            symbols.append(a[i])
        clock_phase += cycles_per_sample
    debug = True
    if debug:
        print("peak frequency index: %d / %d" % (p, len(f)))
        print("samples per symbol: %f" % (1.0/cycles_per_sample))
        print("clock cycles per sample: %f" % (cycles_per_sample))
        print("clock phase in cycles between 1st and 2nd samples: %f" % (clock_phase))
        print("clock phase in cycles at 1st sample: %f" % (clock_phase - cycles_per_sample/2))
        print("symbol count: %d" % (len(symbols)))
    return symbols

def main(Fs, #sampling rate of sdr
        frf, #sdr frequency
        num_samples, #number of samples to get
        debug, #enable/disable debug print/graph-viewing statements
        bps, #expected fsk bps
        fromIQ, #getting samples from iq file or sdr
        IQpath, #path to IQ file
        IQindex #index in IQ file to get samples 
        ):

    if fromIQ:
        f = np.fromfile(open(IQpath), dtype=np.complex64)

        sm = f[IQindex]    

    #get samples from SDR 
    else:
        sdr = SoapySDRDriver.init_radio(device = 'rtlsdr', Fs = 250000, frf = 915e6, gain = None)
        rxStream, buff = SoapySDRDriver.start_rf_stream(sdr, bufsize = num_samples)

        for i in range(5):
            SoapySDRDriver.get_radio_samples(sdr, rxStream, buff, printstats = False)

        #sdr = init_radio(Fs, frf, gain = 5)
        #call this function but don't recieve any samples to get sdr into steady state
        #get_radio_samples(sdr, num_samples)

    #maybe replace with a better fm detection algorithm 
        '''
        signal = np.array([])
        while signal.size == 0:
            sm = get_radio_samples(sdr, num_samples) 
            Fsspect = Fs/30
            num_taps = 51
            h = scipy.signal.firwin(num_taps, Fsspect, nyq=Fs/2)
            smfilt = scipy.signal.lfilter(h, 1.0, sm)

            signal = smfilt
            #smfilt = sm
            psd(smfilt, NFFT=1024, Fs=sdr.sample_rate/1e6, Fc=sdr.center_freq/1e6)
            show()
            plt.plot(10*np.log(np.abs(np.fft.fftshift(np.fft.fft(sm, 2**10)))))
            plt.plot(10*np.log(np.abs(np.fft.fftshift(np.fft.fft(smfilt, 2**10)))))
            plt.show()
            #signal = detectsignal(smfilt, airtime = .1, Fs = Fs, resolution = 10, fftpow = 10)
        '''

        sm = SoapySDRDriver.get_radio_samples(sdr, rxStream, buff, printstats = False)

        if debug:
            SoapySDRDriver.simple_psd(buff, NFFT = 1024, Fs = 250000, frf = 915e6)

        SoapySDRDriver.stop_rf_stream(sdr, rxStream)

        #sm = get_radio_samples(sdr, num_samples) 

    if debug:
        specto(sm, 6, graph = True, fft_sizecust = 2**7)

    #low pass filter to isolate baseband fsk signal
    Fsspect = Fs/3
    num_taps = 51
    h = scipy.signal.firwin(num_taps, Fsspect, nyq=Fs/2)
    smfilt = scipy.signal.lfilter(h, 1.0, sm)

    if debug:
        psd(smfilt, NFFT=1024, Fs=Fs/1e6, Fc=frf/1e6)
        show()

    if debug:
        freqsignal = frequency_detector(smfilt, Fs, averaging = 2, graph = True)
    else:
        freqsignal = frequency_detector(smfilt, Fs, averaging = 2, graph = None)

    if debug:
        freqold1 = freqsignal

    #low pass filter on instantaneous frequency signal to smooth out signal
    Fsspect = Fs/8
    num_taps = 51
    h = scipy.signal.firwin(num_taps, Fsspect, nyq=Fs/2)
    freqsignal = scipy.signal.lfilter(h, 1.0, freqsignal)

    if debug:
        freqold2 = freqsignal

    #getting signal mean for thresholding
    freqmean = get_avg(freqsignal, 1000)

    if debug:
        print(freqmean)

    #hysteresis thresholding
    #if data is in threshold it will be zero 
    '''
    tol = 100
    freqsignalclockrecovery = np.zeros_like(freqsignal)
    freqsignalclockrecovery[ freqsignal < freqmean-tol] = 0
    freqsignalclockrecovery[ freqsignal > freqmean+tol] = 1

    freqsignal = np.zeros_like(freqsignal)
    freqsignal[ freqsignal < freqmean-tol] = -5000
    freqsignal[ freqsignal > freqmean+tol] = 5000
    '''

    #thresholding of signal for clock recovery
    freqsignalclockrecovery = np.where(freqsignal > freqmean, 1, 0)

    #thresholding of signal for bps calculation
    freqsignalthreshold = np.where(freqsignal > freqmean, 5000, -5000)

    if debug:
        print(wpcr(freqsignalclockrecovery))

    if debug:
        plt.plot(freqsignalclockrecovery)
        plt.show()

    #getting pulses at rising/falling edge
    freqdelay = IQ_delay(freqsignalthreshold, 1, real = True)
    freqxor = np.bitwise_xor(freqsignalthreshold.astype(int), freqdelay.astype(int)).astype(float)

    if debug:
        print('Number of Samples: ' + str(freqxor.size))

    #converting back to float
    freqxor = np.bitwise_not(freqxor.astype(int)).astype(float)

    if debug:
        plt.plot(freqold2)
        plt.plot(freqold1)
        plt.plot(freqsignal)
        plt.plot(freqxor * 2500)
        plt.show()
        plt.plot(np.abs(np.fft.fftshift(np.fft.fft(freqxor))))
        plt.show()

    #bandpass filter around expected bps tone
    num_taps = 51
    passfreq = [bps - 1000, bps + 1000]
    h = scipy.signal.firwin(num_taps, passfreq, fs = Fs, pass_zero = False, window = 'hamming')
    freqfilt = scipy.signal.lfilter(h, 1.0, freqxor)

    #high pass filter to reduce dc tone amplitude
    num_taps = 51
    h = scipy.signal.firwin(num_taps, cutoff = bps, fs = Fs, pass_zero = 'highpass', window = 'hamming')
    freqfilt = scipy.signal.lfilter(h, 1.0, freqfilt)

    freqfft = np.abs(np.fft.fftshift(np.fft.fft(freqfilt)))

    #position of 9600 Hz peak will be (samplerate/2)*tonelocation/(samplenum/2)
    if debug:
        print('Expected Bin of Bit Rate: ' + str((bps*((freqfft.size/2))/(Fs/2))))
        plt.plot(freqfft[int(freqfft.size/2)::])
        plt.show()

    freqhalf = freqfft[int(freqfft.size/2)::]
    fftargmax = np.argmax(freqhalf)
    #if dc componet is higher than first harmonic delete it and redo argmax
    if fftargmax == 0:
        freqhalf = np.delete(freqhalf, 0)
        fftargmax = np.argmax(freqhalf)

    if debug:
        print('Obtained Bit Rate Bin: ' + str(fftargmax))


    #getting bit rate from fft bin: bps = tonebinnum * (Fs / 2) / (fftsize / 2)
    fskbitrate = (fftargmax) * Fs / freqfft.size
    print('Bit Rate: ' + str(fskbitrate))

    return fskbitrate

#config options

#sampling rate of sdr
#Fs = 250000
#sdr frequency
#frf = 915e6
#number of samples to get
#num_samples = 10000
#enable/disable debug print/graph-viewing statements
#debug =  False
#expected fsk bps
#bps = 4800
#getting samples from iq file or sdr
#fromIQ = True
#path to IQ file
#IQpath = r"C:\Users\jstephenson\Documents\IQ RECORDINGS\LFMFSK4800", 
#indexes in IQ file to process
#IQindex = list(range(390000,420000))



if __name__ == "__main__":
    main(Fs = 250000, 
    frf = 915e6, 
    num_samples = 10000, 
    debug = False, 
    bps = 4800, 
    fromIQ = True, 
    IQpath = r"C:\Users\jstephenson\Documents\IQ RECORDINGS\LFMFSK4800", 
    IQindex = list(range(390000,420000)))