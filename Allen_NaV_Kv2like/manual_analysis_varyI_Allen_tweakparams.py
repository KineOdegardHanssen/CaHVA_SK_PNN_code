import numpy
import matplotlib.pyplot as plt

def avg_and_rms(x):
    N = len(x)
    avgx = numpy.mean(x)
    rmsx = 0
    for i in range(N):
        rmsx += (avgx-x[i])**2
    rmsx = numpy.sqrt(rmsx/(N-1))
    return avgx,rmsx

def manual(filename,idelay,idur,spikedurat,skiptime):
    """Manual approach"""

    # Use numpy to read the trace data from the txt file
    data = numpy.loadtxt(filename)

    # Time is the first column
    time = data[:, 0]
    # Voltage is the second column
    voltage = data[:, 1]
    
    anfraction = (idur-skiptime)/float(idur)
    tstartan   = idelay+skiptime
    print('anfraction:',anfraction)
    print('tstartan:',tstartan)
    vmax = max(voltage) 
    vmin = min(voltage) 
    deltav = vmax-vmin
    vthr  = -20   # If there is a peak above this value, we count it 
    vprev = vthr-40 # A peak never kicks in at initiation, even if I change vthr
    durthr = spikedurat # Height at which we measure the duration. 
    Npeaks = 0
    peakmins  = []
    peakvals  = []
    peaktimes = []
    passtimes_up = []
    passvals_up  = []
    passtimes_down = []
    passvals_down  = []
    minalready = False
    for i in range (1,len(voltage)-1):  
        if time[i]<idelay+idur:
            if voltage[i-1]<voltage[i] and voltage[i+1]<voltage[i] and voltage[i]>vthr and time[i]>tstartan:
                peaktimes.append(time[i])
                peakvals.append(voltage[i])
                Npeaks+=1
                minalready = False
            if voltage[i-1]>voltage[i] and voltage[i+1]>voltage[i] and voltage[i]<vthr and time[i]>tstartan and minalready==False and Npeaks>0:
                peakmins.append(voltage[i])
                minalready = True
            if voltage[i]>=durthr and voltage[i-1]<durthr and time[i]>tstartan: # Passing upwards
                tbef = time[i-1]
                taft = time[i]
                Vbef = voltage[i-1]
                Vaft = voltage[i]
                a = (Vaft-Vbef)/(taft-tbef)
                b = Vbef-a*tbef
                tint = (durthr-b)/a
                Vint = a*tint+b
                passtimes_up.append(tint)
                passvals_up.append(Vint) # For plotting
            elif voltage[i]>=durthr and voltage[i+1]<durthr and len(passtimes_up)>0: # Passing downwards
                tbef = time[i]
                taft = time[i+1]
                Vbef = voltage[i]
                Vaft = voltage[i+1]
                a = (Vaft-Vbef)/(taft-tbef)
                b = Vbef-a*tbef
                tint = (durthr-b)/a
                Vint = a*tint+b
                passtimes_down.append(tint)
                passvals_down.append(Vint) # For plotting
        else:
            break
    
    Npeaks /= anfraction # Want frequency
    # Checking if we've got consistent firing:
    if Npeaks!=0:
        if peaktimes[-1]<=(idur/2.+idelay): #Checking if there's no firing in the last half of the stim. interval
            Npeaks=0                        
    
    dur = []
    isi = []
    amps = []
    Namps = min([len(peakmins),len(peakvals)])
    Ndur = min([len(passtimes_up),len(passtimes_down)]) 
    for i in range(Ndur-1):
        dur.append(passtimes_down[i]-passtimes_up[i])
    for i in range(Namps):
        amps.append(peakvals[i]-peakmins[i])
    for i in range(1,len(peaktimes)):
        isi.append(peaktimes[i]-peaktimes[i-1])
    time_peakvals = peaktimes
    
    '''
    plt.figure(figsize=(6,5))
    plt.plot(time,voltage,',')
    plt.plot(time_peakvals,peakvals,'o',label='peaks')
    plt.plot(passtimes_up,passvals_up,'o',label='dur basis, up')
    plt.plot(passtimes_down,passvals_down,'o',label='dur basis, down')
    plt.xlabel('Time [ms]')
    plt.ylabel('Voltage [mV]')
    plt.title('Testing implementation')
    plt.legend(loc='lower left')
    plt.tight_layout()
    plt.show() 
    '''
    
    ## Avg and rms:
    amps_avg, amps_rms = avg_and_rms(amps)
    peakmins_avg, peakmins_rms = avg_and_rms(peakmins)
    peakvals_avg, peakvals_rms = avg_and_rms(peakvals)
    dur_avg, dur_rms = avg_and_rms(dur)
    isi_avg, isi_rms = avg_and_rms(isi)
    
    return Npeaks, peaktimes, peakmins_avg, peakmins_rms, peakvals_avg,  peakvals_rms, dur_avg, dur_rms, isi_avg, isi_rms, isi, dur, amps_avg, amps_rms

if __name__ == '__main__':
    skiptime   = 500 
    spikedurat = -40
    idur       = 1000 # ms
    idelay     = 100
    v_init     = -86 # mV
    Ra         = 100
    somasize   = 10  
    dtexp      = -7
    t_before_rec = -600.
    
    model_folder = '' # We are in the model folder
    
    cm = 1.0
    iamps = [0,0.001,0.002,0.003,0.004,0.005,0.006,0.007,0.008,0.009,0.01,0.011,0.012,0.013,0.014,0.015,0.016,0.017,0.018,0.019,0.02,0.021,0.022,0.023,0.024,0.025,0.026,0.027,0.028,0.029,0.03,0.031,0.032,0.033,0.034,0.035,0.036,0.037,0.038,0.039,0.04,0.041,0.042,0.043,0.044,0.045,0.046,0.047,0.048,0.049]#[0,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.11,0.12,0.13,0.14,0.15,0.16,0.17,0.18,0.19,0.2,0.21,0.22,0.23,0.24,0.25,0.26,0.27,0.28,0.29,0.3,0.31,0.32,0.33,0.34,0.35,0.36,0.37,0.38,0.39,0.4,0.41,0.42,0.43,0.44,0.45,0.46,0.47,0.48,0.49,0.5,0.51,0.52,0.53,0.54,0.55,0.56,0.57,0.58,0.59,0.6,0.61,0.62,0.63,0.64,0.65,0.66,0.67,0.68,0.69,0.7,0.71,0.72,0.73,0.74,0.75,0.76,0.77,0.78,0.79,0.8]
    Namps = len(iamps)
    
    varyE        = 0     
    varymech     = 'Epas'
    vary_NaV     = True
    vary_Kv2like = True
    
    gnav     = 0.5
    gkv2like = 0.5

    namestring = ''
    if varymech=='ENa':
        varyE = 63 
        namestring = namestring + 'ENa'+str(varyE)
    elif varymech=='EK':
        varyE = -97
        namestring = namestring + 'EK'+str(varyE)
    elif varymech=='Epas': 
        varyE = -20 # Vary by shifts
        namestring = namestring + 'Epasshift'+str(varyE)
    if vary_NaV==True:
        namestring = namestring + '_gNaV'+str(gnav)+'p'
    if vary_Kv2like==True:
        namestring = namestring + '_gKv2like'+str(gkv2like)+'p'
    namestring = namestring +'_'
    
    
    Nspikes = numpy.zeros(Namps)
    avg_ISI = numpy.zeros(Namps)
    rms_ISI = numpy.zeros(Namps)
    avg_AP_ampl = numpy.zeros(Namps)
    rms_AP_ampl = numpy.zeros(Namps)
    avg_AP_mins = numpy.zeros(Namps)
    rms_AP_mins = numpy.zeros(Namps)
    avg_AP_halfwidth = numpy.zeros(Namps)
    rms_AP_halfwidth = numpy.zeros(Namps)
    avg_AP_amplitudes = numpy.zeros(Namps)
    rms_AP_amplitudes = numpy.zeros(Namps)
    
    # Set names
    outfolder = 'Results/Soma%i/' % somasize
    outfilename_Nspikes = outfolder+'somaAllen_idur%i_varyiamp'% (idur)+'_manual_cm'+str(cm)+'_'+namestring+'Nspikes_vs_I_s'+str(skiptime)+'.txt'
    outfilename_amps    = outfolder+'somaAllen_idur%i_varyiamp'% (idur)+'_manual_cm'+str(cm)+'_'+namestring+'Ampl_vs_I_s'+str(skiptime)+'.txt'
    outfilename_APampl  = outfolder+'somaAllen_idur%i_varyiamp'% (idur)+'_manual_cm'+str(cm)+'_'+namestring+'Vmax_vs_I_s'+str(skiptime)+'.txt'
    outfilename_APmins  = outfolder+'somaAllen_idur%i_varyiamp'% (idur)+'_manual_cm'+str(cm)+'_'+namestring+'Vmin_vs_I_s'+str(skiptime)+'.txt'
    outfilename_APdhw   = outfolder+'somaAllen_idur%i_varyiamp'% (idur)+'_manual_cm'+str(cm)+'_'+namestring+'sdurat%s' % str(spikedurat)+'_vs_I_s'+str(skiptime)+'.txt'
    outfilename_ISI     = outfolder+'somaAllen_idur%i_varyiamp'% (idur)+'_manual_cm'+str(cm)+'_'+namestring+'ISI_vs_I_s'+str(skiptime)+'.txt'
    # make files
    outfile_Nspikes = open(outfilename_Nspikes,'w')
    outfile_amps    = open(outfilename_amps,'w')
    outfile_APampl  = open(outfilename_APampl,'w')
    outfile_APmins  = open(outfilename_APmins,'w')
    outfile_APdhw   = open(outfilename_APdhw,'w')
    outfile_ISI     = open(outfilename_ISI,'w')
    firedbefore = False
    for j in range(Namps):
        iamp = iamps[j]
        print('Step ', j+1, ' of', Namps, 'cm:', cm, 'iamp:', iamp)
        infolder = outfolder + 'current_idur%i_iamp'%idur+str(iamp)+'/'
        filename = infolder+namestring+'somaonly_cm'+str(cm)+'_idur%i_iamp'%idur+str(iamp)+'_dtexp%i_V.txt' % dtexp
        Nspikes[j], peaktimes, avg_AP_mins[j], rms_AP_mins[j], avg_AP_ampl[j], rms_AP_ampl[j], avg_AP_halfwidth[j], rms_AP_halfwidth[j], avg_ISI[j], rms_ISI[j], ISI, durs, avg_AP_amplitudes[j], rms_AP_amplitudes[j] = manual(filename,idelay,idur,spikedurat,skiptime)
        if Nspikes[j]!=0 and firedbefore==False:
            firedbefore = True    
        if Nspikes[j]!=0 or firedbefore==False:
            outfile_Nspikes.write('%.5f %i\n' % (iamp,Nspikes[j]))
        if Nspikes[j]!=0:
            outfile_amps.write('%.5f %.10f %.10f\n' % (iamp,avg_AP_amplitudes[j],rms_AP_amplitudes[j]))
            outfile_APampl.write('%.5f %.10f %.10f\n' % (iamp,avg_AP_ampl[j],rms_AP_ampl[j]))
            outfile_APmins.write('%.5f %.10f %.10f\n' % (iamp,avg_AP_mins[j],rms_AP_mins[j]))
            outfile_APdhw.write('%.5f %.10f %.10f\n' % (iamp,avg_AP_halfwidth[j],rms_AP_halfwidth[j]))
            # Write all durs:
            outfilename_durs_all = infolder+'somaAllen_idur%i_iamp'% (idur)+str(iamp) +'_manual_cmf'+str(cm)+'_alldurs_s'+str(skiptime)+'.txt'
            figname_durs_all     = infolder+'somaAllen_idur%i_iamp'% (idur)+str(iamp) +'_manual_cmf'+str(cm)+'_alldurs'+str(skiptime)+'.png'
            outfile_durs_all = open(outfilename_durs_all,'w')
            for k in range(len(durs)):    
                outfile_durs_all.write('%.10f ' % durs[k])
            outfile_durs_all.close()
            ptimes = peaktimes[1:-1]
            if len(ptimes)==len(durs):
                plt.figure(figsize=(6,5))
                plt.plot(ptimes,durs)
                plt.xlabel('Time [ms]')
                plt.ylabel('Voltage [mV]')    
                plt.title('I=%s nA' % str(iamp))
                plt.tight_layout()
                plt.savefig(figname_durs_all)
            # Write all peak times (matching durs):    
            outfilename_pt_all = infolder+'somaAllen_idur%i_iamp'% (idur)+str(iamp) +'_manual_cmf'+str(cm)+str(v_init)+'_trec'+str(t_before_rec)+'_'+namestring+'_peaktimes'+str(skiptime)+'.txt'
            outfile_pt_all = open(outfilename_pt_all,'w')
            for k in range(len(peaktimes)-1):    
                outfile_pt_all.write('%.10f ' % peaktimes[k])
            outfile_pt_all.close()
            outfile_ISI.write('%.5f %.10f %.10f\n' % (iamp,avg_ISI[j],rms_ISI[j]))
            # Write all ISIs:
            outfilename_ISI_all = infolder+'somaAllen_idur%i_iamp'% (idur)+str(iamp) +'_manual_cmf'+str(cm)+str(v_init)+'_trec'+str(t_before_rec)+'_'+namestring+'ISIall'+str(skiptime)+'.txt'
            outfile_ISI_all = open(outfilename_ISI_all,'w')
            for k in range(len(ISI)):
                outfile_ISI_all.write('%.10f ' % ISI[k])
            outfile_ISI_all.close()
    outfile_Nspikes.close()
    outfile_amps.close()
    outfile_APampl.close()
    outfile_APmins.close()
    outfile_APdhw.close()
    outfile_ISI.close()

