import os
from os.path import join
import sys
import matplotlib.pyplot as plt
import json
import neuron
import time as tm
import numpy as np

import LFPy

def return_allen_cell_model(model_folder, cm_factor, somasize, idur,iamp,idelay,tstop,dtexp,gnaf,gkaf,gcahva,gpas,gbk,vary_naf,vary_kaf,vary_Ca_HVA,vary_bk,vary_gpas,varyE,varymech,v_init,t_before_rec):

    params = json.load(open(join(model_folder, "fit_parameters.json"), 'r'))

    Ra = params["passive"][0]["ra"]
    
    e_pas = params["passive"][0]["e_pas"]
    celsius = params["conditions"][0]["celsius"]
    cm = params["passive"][0]["cm"]
    reversal_potentials = params["conditions"][0]["erev"]
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
        'v_init': v_init,    # initial membrane potential 
        'passive': False,   # turn on NEURONs passive mechanism for all sections
        'nsegs_method': 'fixed_length', # spatial discretization method
        'max_nsegs_length': 20.,
        'dt': mydt,      # simulation time step size
        'tstart': t_before_rec,    # start time of simulation, recorders start at t=0 #
        'tstop': tstop,     # stop simulation at this time.
    }
    
    cell = LFPy.Cell(**cell_params)
    
    
    for sec in neuron.h.allsec():
        sec.insert("pas")
        sectype = sec.name().split("[")[0]
        for sec_dict in active_mechs:
            if sec_dict["section"] == sectype:
                if sectype=="soma":
                    sec.cm   = cm*cm_factor
                    sec.L    = somasize
                    sec.diam = somasize
                if not sec_dict["mechanism"] == "":
                    print('sec_dict[mechanism]:',sec_dict["mechanism"])
                    sec.insert(sec_dict["mechanism"])
                exec("sec.{} = {}".format(sec_dict["name"], sec_dict["value"]))
    
    for sec in neuron.h.allsec():
        sectype = sec.name().split("[")[0]
        # First: Change mechanisms 
        # Mechanisms that can be anywhere
        if varymech=='Epas':
            sec.e_pas += varyE
        if sectype=='soma':
            if varymech=='ENa':
                sec.ena = varyE
            elif varymech=='EK':
                sec.ek = varyE
            if vary_naf==True:
                sec.gbar_naf *= gnaf
            if vary_kaf==True:
                sec.gbar_kaf *= gkaf
            if vary_bk==True:
                sec.gkbar_SK2 *= gbk
            if vary_Ca_HVA==True:
                sec.pbar_caq *= gcahva
            if vary_gpas==True: 
                sec.g_pas *= gpas
    
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
    iamps =  [0,0.001,0.002,0.003,0.004,0.005,0.006,0.007,0.008,0.009,0.01,0.011,0.012,0.013,0.014,0.015,0.016,0.017,0.018,0.019,0.02,0.021,0.022,0.023,0.024,0.025,0.026,0.027,0.028,0.029,0.03,0.031,0.032,0.033,0.034,0.035,0.036,0.037,0.038,0.039,0.04,0.041,0.042,0.043,0.044,0.045,0.046,0.047,0.048,0.049,0.05,0.051,0.052,0.053,0.054,0.055,0.056,0.057,0.058,0.059,0.06]
    idur   = 1000 
    idelay = 100
    tstop  = idur+idelay+10.
    dtexp  = -7
    v_init = -86
    t_before_rec = -600.
    Ncm = len(cm_factors)
    Ni  = len(iamps)
    
    somasize = 10
    
    model_folder = '' # We are in the model folder
    
    
    vary_naf    = False
    vary_kaf    = False
    vary_Ca_HVA = True # False
    vary_bk     = True # False
    vary_gpas   = False #
    
    gnaf    = 1.0
    gkaf    = 1.0
    gcahvas = [0.2]
    gbks    = [15]
    gpas    = 1.0
    
    for gcahva in gcahvas:
        for gbk in gbks:
            varyE = 0
            varymech = 'None' # 'Epas' # 'EK' # 'ENa' # 
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
        
            if vary_naf==True:
                namestring = namestring + '_gnaf'+str(gnaf)+'p'
            if vary_kaf==True:
                namestring = namestring + '_gkaf'+str(gkaf)+'p'
            if vary_bk==True:
                namestring = namestring + '_gbk'+str(gbk)+'p'
            if vary_Ca_HVA==True:
                namestring = namestring + '_gCaHVA'+str(gcahva)+'p'
            if vary_gpas==True: 
                namestring = namestring + '_gpas'+str(gpas)+'p'
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
                    t, v = return_allen_cell_model(model_folder,cm_factor,somasize,idur,iamp,idelay,tstop,dtexp,gnaf,gkaf,gcahva,gpas,gbk,vary_naf,vary_kaf,vary_Ca_HVA,vary_bk,vary_gpas,varyE,varymech,v_init,t_before_rec)
                    
                    # Save results:
                    filename = folder+namestring+'somaonly_cm'+str(cm_factor)+'_idur%i_iamp'%idur+str(iamp)+'_dtexp%i_vinit' % dtexp+str(v_init)+'_trec'+str(t_before_rec)+'_V.txt'
                    plotname   = folder+namestring+'somaonly_cm'+str(cm_factor)+'_idur%i_iamp'%idur+str(iamp)+'_dtexp%i_vinit' % dtexp+str(v_init)+'_trec'+str(t_before_rec)+'V.png'    
                    
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
                    #plt.show()
                    
                    print('v[-1]:',v[-1])
                    
