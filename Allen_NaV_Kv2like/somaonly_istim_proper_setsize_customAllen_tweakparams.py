import os
from os.path import join
import sys
import matplotlib.pyplot as plt
import json
import neuron
import time as tm
import numpy as np

import LFPy

def return_allen_cell_model(model_folder, cm_factor, somasize, idur,iamp,idelay,tstop,dtexp,gnav,gkv2like,vary_NaV,vary_Kv2like,varymech,varyE):

    params = json.load(open(join(model_folder, "fit_parameters.json"), 'r'))

    Ra = params["passive"][0]["ra"]
    
    e_pas = params["passive"][0]["e_pas"]
    celsius = params["conditions"][0]["celsius"]
    cm = params["passive"][0]["cm"]
    reversal_potentials = params["conditions"][0]["erev"]
    v_init = params["conditions"][0]["v_init"]
    active_mechs = params["genome"]
    Ra = params["passive"][0]["ra"] 
    neuron.h.celsius = celsius
    cm = cm[0]["cm"]
    
    print('cm:',cm)
    print('e_pas:',e_pas)

    mydt = 2.**dtexp
    # Define cell parameters
    cell_params = {
        'morphology': "somaonly.hoc",
        'v_init': -86.8,    # initial membrane potential
        'passive': False,   # turn on NEURONs passive mechanism for all sections
        'nsegs_method': 'fixed_length', # spatial discretization method
        'max_nsegs_length': 20.,
        'dt': mydt,      # simulation time step size
        'tstart': -600.,    # start time of simulation, recorders start at t=0
        'tstop': tstop,     # stop simulation at this time.
    }
    
    cell = LFPy.Cell(**cell_params)
    
    
    for sec in neuron.h.allsec():
        sec.insert("pas")
        sectype = sec.name().split("[")[0]
        for sec_dict in active_mechs:
            if sec_dict["section"] == sectype:
                if sectype=="soma": # Works
                    sec.cm   = cm*cm_factor
                    sec.L    = somasize
                    sec.diam = somasize
                if not sec_dict["mechanism"] == "":
                    sec.insert(sec_dict["mechanism"])
                exec("sec.{} = {}".format(sec_dict["name"], sec_dict["value"]))
    
    for sec in neuron.h.allsec():
        sectype = sec.name().split("[")[0]
        # First: Change mechanisms 
        # Mechanisms that can be anywhere
        if varymech=='pas':
            if varyE!='None':
                sec.e_pas += varyE
        if sectype=='soma':
            if varymech=='Na':
                if varyE!='None':
                    sec.ena = varyE
            elif varymech=='K':
                if varyE!='None':
                    sec.ek = varyE
            if vary_NaV==True:
                sec.gbar_NaV *= gnav
            if vary_Kv2like==True:
                sec.gbar_Kv2like *= gkv2like
    
    stim_idx = 0
    stim_params = {
                'idx': stim_idx,
                'record_current': True,
                'pptype': 'IClamp',
                'amp': iamp,
                'dur': idur,
                'delay': idelay,
            }
    
    stimulus = LFPy.StimIntElectrode(cell, **stim_params)
    #syn.set_spike_times(np.array([1]))
    cell.simulate(rec_vmem=True, rec_imem=True)
    t, v = cell.tvec.copy(), cell.vmem[0].copy()
    
    if cm_factor==1 and iamp==0.02:
        for sec in neuron.h.allsec():
            neuron.h.psection()
    
    cell.__del__()
    return t, v


if __name__ == '__main__':
    cm_factors = [1.0]
    iamps = [0.01]
    idur   = 100
    idelay = 10
    tstop  = idur+idelay+20.
    dtexp  = -7
    Ncm = len(cm_factors)
    Ni  = len(iamps)
    
    somasize = 10
    
    model_folder = '' # We are in the model folder
    
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
    
    print('iamps:',iamps)
    for i in range(Ni):
        iamp = iamps[i]
        folder = 'Results/Soma%i/current_idur%i_iamp'%(somasize,idur)+str(iamp)+'/'
        folder_total = folder
        if os.path.exists(folder_total)==False:
            os.mkdir(folder_total)
        print('Step ', i+1, ' out of ', Ni)
        for j in range(Ncm):
            cm_factor = cm_factors[j]
            t, v = return_allen_cell_model(model_folder,cm_factor,somasize,idur,iamp,idelay,tstop,dtexp,gnav,gkv2like,vary_NaV,vary_Kv2like,varymech,varyE)
            
            # Save results:
            filename = folder+namestring+'somaonly_cm'+str(cm_factor)+'_idur%i_iamp'%idur+str(iamp)+'_dtexp%i_V.txt' % dtexp
            plotname   = folder+namestring+'somaonly_cm'+str(cm_factor)+'_idur%i_iamp'%idur+str(iamp)+'_dtexp%i_V.png' % dtexp
            
            Nt = len(t)
            file = open(filename,'w')
            for k in range(Nt):
                file.write('%.16e %.16e\n' % (t[k],v[k]))
            file.close()   
            
            plt.figure(figsize=(6,5))
            plt.plot(t,v)
            plt.xlabel('Time (ms)')
            plt.ylabel('Membrane potential (mV)')
            plt.title('Voltage vs time, %.2f*Cm' % cm_factor)
            plt.savefig(plotname)
    
