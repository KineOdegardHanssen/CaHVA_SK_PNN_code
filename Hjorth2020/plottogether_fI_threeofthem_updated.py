import numpy
import matplotlib.pyplot as plt
from matplotlib import gridspec


def avg_and_rms(x):
    N = len(x)
    avgx = numpy.mean(x)
    rmsx = 0
    for i in range(N):
        rmsx += (avgx-x[i])**2
    rmsx = numpy.sqrt(rmsx/(N-1))
    return avgx,rmsx

vshift      = 0
spikedurat = -40
idur       = 1000 # ms
idelay     = 10
v_init     = -86 # mV
Ra         = 100
somasize   = 10
dtexp      = -7
skiptime   = 500
t_before_rec = -600.

cm = 1.0

showplot = True # False

plt.rc('xtick',labelsize=14)
plt.rc('ytick',labelsize=14)
    
fig = plt.figure(figsize=(15, 5),dpi=300)
    
gs = gridspec.GridSpec(2, 6, height_ratios=[4, 1])
    
ax1 = plt.subplot(gs[0, 0:2])
ax2 = plt.subplot(gs[0, 2:4])
ax3 = plt.subplot(gs[0, 4:6])
ax4 = plt.subplot(gs[1, 0:6])

ax1.set_title(r'A',loc='left',fontsize=18)
ax2.set_title(r'B',loc='left',fontsize=18)
ax3.set_title(r'B',loc='left',fontsize=18)
ax1.set_title(r'Altering $\bar{g}_\mathregular{CaHVA}$, with $0\bar{g}_\mathregular{SK}$',fontsize=18)
ax2.set_title(r'Altering $\bar{g}_\mathregular{SK}$',fontsize=18)
ax3.set_title(r'Altering $\bar{g}_\mathregular{CaHVA}$, with 1.0$\bar{g}_\mathregular{SK}$',fontsize=18)
    
ax1.set_xlabel(r'$I$ [nA]',fontsize=14)
ax1.set_ylabel(r'$f$ [Hz]',fontsize=14)
ax2.set_xlabel(r'$I$ [nA]',fontsize=14)
ax2.set_ylabel(r'$f$ [Hz]',fontsize=14)
ax3.set_xlabel(r'$I$ [nA]',fontsize=14)
ax3.set_ylabel(r'$f$ [Hz]',fontsize=14)

varyE = 0
varymech = 'None' # 'Epas' # 'EK' # 'ENa' # 
outnamestring = 'vary'
if varymech=='ENa':
    varyE = 63
    outnamestring = outnamestring + '_ENa'+str(varyE)
elif varymech=='EK':
    varyE = -97#-107
    outnamestring = outnamestring + '_EK'+str(varyE)
elif varymech=='Epas':
    varyE = -20 # Vary by shifts
    outnamestring = outnamestring + '_Epasshift'+str(varyE)


vary_naf     = False
vary_kaf     = False
vary_Ca_HVA  = True # False
vary_gpas    = False # 

gnafs   = 1.0
gkafs   = 1.0
gcahvas = [0.02,0.04,0.06,0.08,0.1,0.12,0.14,0.16,0.18,0.2,0.25,0.3,0.4,0.5,1.0,2.0,2.5]
gcahvas_plot = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.25,1.5,2.0,2.5,5.0,10.0,12.5]
gpases  = 1.0
Ng      = len(gcahvas)

colors = []
for i in range(Ng):
    colors.append((1.0-i/float(Ng),0,i/float(Ng),1.0-abs(0.5-i/float(Ng))))

if vary_naf==True:
    outnamestring = outnamestring + '_gnaf'
if vary_kaf==True:
    outnamestring = outnamestring + '_gkaf'
if vary_Ca_HVA==True:
    outnamestring = outnamestring + '_gCaHVA'
if vary_gpas==True: 
    outnamestring = outnamestring + '_gpas'
outnamestring = outnamestring +'_'

folder    = 'CaHVA_Allen/Shift_gCaHVA_V/Vshift_'+str(vshift)+'/Results/Soma%i/' % somasize
outfolder = 'Compare/Soma10/'
figname   = outfolder+'mainmodel_threeofI.png'

i = 0
for gcahva in gcahvas:
    gcahva_plot = gcahvas_plot[i]
    namestring = ''
    thelabel = ''
    if vary_naf==True:
        namestring = namestring + '_gnaf'+str(gnaf)+'p'
    if vary_kaf==True:
        namestring = namestring + '_gkaf'+str(gkaf)+'p'
    if vary_Ca_HVA==True:
        namestring = namestring + '_gCaHVA'+str(gcahva)+'p'    
        thelabel   = thelabel + r'%s$\bar{g}$' % str(gcahvas_plot[i])
    if vary_gpas==True: 
        namestring = namestring + '_gpas'+str(gpas)+'p'
    namestring = namestring +'_'
    
    # Set names
    infilename_Nspikes = folder+'somaHjorth_idur%i_varyiamp'% (idur)+'_manual_cm'+str(cm)+namestring+'_Nspikes_vs_I_s'+str(skiptime)+'.txt'
    # make files    
    infile_Nspikes = open(infilename_Nspikes,'r')
    lines_Nspikes = infile_Nspikes.readlines()
    
    I_Nspikes = []
    Nspikes   = []
    
    for line in lines_Nspikes:
        words = line.split()
        if len(words)>0:
            I_Nspikes.append(float(words[0]))
            Nspikes.append(float(words[1]))
    
    infile_Nspikes.close()
    
    if gcahva_plot==1.0:
        thislinewidth = 3
    else:
        thislinewidth = 2
    
    if numpy.sum(Nspikes)!=0:
        ax1.plot(I_Nspikes,Nspikes, label=thelabel,linewidth=thislinewidth,color=colors[i])
    i+=1


vary_SK      = True # False # 
vary_Ca_HVA  = True # False
    
gnafs   = 1.0
gkafs   = 1.0
gsks    = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.25,1.5,2.0,2.5,5.0,10.0,12.5]
gcahva  = 0.2
gpases  = 1.0
Ng      = len(gsks)
    
colors = []
for i in range(Ng):
    colors.append((1.0-i/float(Ng),0,i/float(Ng),1.0-abs(0.5-i/float(Ng))))

folder    = 'CaHVA_Allen/SK_Allen/Shift_gCaHVA_V/Vshift_'+str(vshift)+'/Results/Soma%i/' % somasize    
i = 0
for gsk in gsks:
    namestring = ''
    thelabel = ''
    if vary_naf==True:
        namestring = namestring + '_gnaf'+str(gnaf)+'p'
    if vary_kaf==True:
        namestring = namestring + '_gkaf'+str(gkaf)+'p'  
    if vary_SK==True:
        namestring = namestring + '_gSK'+str(gsk)+'p'
        thelabel   = thelabel + r'%s$\bar{g}$' % str(gsk)
    if vary_Ca_HVA==True:
        namestring = namestring + '_gCaHVA'+str(gcahva)+'p'
    if vary_gpas==True: 
        namestring = namestring + '_gpas'+str(gpas)+'p'
    namestring = namestring +'_'
    
    # Set names
    infilename_Nspikes = folder+'somaHjorth_idur%i_varyiamp'% (idur)+'_manual_cm'+str(cm)+namestring+'_Nspikes_vs_I_s'+str(skiptime)+'.txt'
    # make files    
    infile_Nspikes = open(infilename_Nspikes,'r')
    lines_Nspikes = infile_Nspikes.readlines()
    
    I_Nspikes = []
    Nspikes   = []
    
    for line in lines_Nspikes:
        words = line.split()
        if len(words)>0:
            I_Nspikes.append(float(words[0]))
            Nspikes.append(float(words[1]))
    
    infile_Nspikes.close()
    
    if gsk==1.0:
        thislinewidth = 3
    else:
        thislinewidth = 2
    
    if numpy.sum(Nspikes)!=0:
        ax2.plot(I_Nspikes,Nspikes, label=thelabel,linewidth=thislinewidth,color=colors[i])
    i+=1

########################################################################
vary_SK      = True # False # 
vary_Ca_HVA  = True # False
    
gnafs   = 1.0
gkafs   = 1.0
gsk     = 1.0
gcahvas = [0.02,0.04,0.06,0.08,0.1,0.12,0.14,0.16,0.18,0.2,0.25,0.3,0.4,0.5,1.0,2.0,2.5]
gcahvas_plot = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.25,1.5,2.0,2.5,5.0,10.0,12.5]
gpases  = 1.0
Ng      = len(gcahvas)
    
colors = []
for i in range(Ng):
    colors.append((1.0-i/float(Ng),0,i/float(Ng),1.0-abs(0.5-i/float(Ng))))

folder    = 'CaHVA_Allen/SK_Allen/Results/Soma%i/' % somasize    
i = 0
for gcahva in gcahvas:
    namestring = ''
    thelabel = ''
    if vary_naf==True:
        namestring = namestring + '_gnaf'+str(gnaf)+'p'
    if vary_kaf==True:
        namestring = namestring + '_gkaf'+str(gkaf)+'p' 
    if vary_SK==True:
        namestring = namestring + '_gSK'+str(gsk)+'p'
    if vary_Ca_HVA==True:
        namestring = namestring + '_gCaHVA'+str(gcahva)+'p'    
        thelabel   = thelabel + r'%s$\bar{g}$' % str(gcahvas_plot[i])
    if vary_gpas==True: 
        namestring = namestring + '_gpas'+str(gpas)+'p'
    namestring = namestring +'_'
    
    # Set names
    infilename_Nspikes = folder+'somaHjorth_idur%i_varyiamp'% (idur)+'_manual_cm'+str(cm)+namestring+'_Nspikes_vs_I_s'+str(skiptime)+'.txt'
    # make files    
    infile_Nspikes = open(infilename_Nspikes,'r')
    lines_Nspikes = infile_Nspikes.readlines()
    
    I_Nspikes = []
    Nspikes   = []
    
    for line in lines_Nspikes:
        words = line.split()
        if len(words)>0:
            I_Nspikes.append(float(words[0]))
            Nspikes.append(float(words[1]))
    
    infile_Nspikes.close()
    
    if gcahva==0.2:
        thislinewidth = 3
    else:
        thislinewidth = 2
    
    if numpy.sum(Nspikes)!=0:
        ax3.plot(I_Nspikes,Nspikes, label=thelabel,linewidth=thislinewidth,color=colors[i])
    i+=1

ax4.set_frame_on(False)
ax4.get_xaxis().set_visible(False)
ax4.get_yaxis().set_visible(False)
for i in range(Ng):
    if gcahvas_plot[i]==1.0:
        thislinewidth = 3
    else:
        thislinewidth = 2
    ax4.plot('-',color=colors[i],linewidth=thislinewidth,label=r'%.1f$\bar{g}$' % gcahvas_plot[i])
ax4.legend(loc='upper center',ncol=9, fontsize=12.5)

plt.tight_layout()
plt.savefig(figname)
if showplot==True:
    plt.show()
