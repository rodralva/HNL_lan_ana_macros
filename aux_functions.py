import uproot
import numpy as np
import awkward as ak

SBND_POTs   = 1e21

def get_true_signal_in_all_spills(ev,pot_scale,DEBUG=False):
    nSig    =sum(ak.flatten(ev["nu_event_type"].array())==0)*pot_scale
    nNonFV  =sum(ak.flatten(ev["nu_event_type"].array())==1)*pot_scale
    
    if DEBUG:
        print("---------------")
        print("True signal in all spills: ", nSig)
        print("True non-FV in all spills: ", nNonFV)
        print("True total in all spills: ", nSig+nNonFV)
    
    return nSig, nNonFV

def get_reco_signal_in_all_spills(ev,pot_scale,DEBUG=False):
    
    ev_type=ak.flatten(ev['slc_true_event_type'].array() )
    ev_comp=ak.flatten(ev['slc_comp']           .array() )>0.5
    #compute comp only once
    
    nSig    =sum((ev_type==0)*ev_comp)*pot_scale
    nNonFV  =sum((ev_type==1)*ev_comp)*pot_scale
    
    if DEBUG:
        print("---------------")
        print("True signal in all spills: ", nSig)
        print("True non-FV in all spills: ", nNonFV)
        print("True total in all spills: ", nSig+nNonFV)
    
    return nSig, nNonFV

event_label_dict = {                
                0:  r"Signal HNL ${\nu e^+e^-}$"       
                , 1:  r"Non-FV HNL ${\nu e^+e^-}$"        
                , 2:  "Dirt HNL"         
                , 3:  "NC ${\pi}^{0}$"        
                , 4:  "Other NC"         
                , 5:  "CC $\\nu_{\mu}$"          
                , 6:  "CC $\\nu_{e}$"           
                , 7:  "Non-FV $\\nu$"         
                , 8:  "Dirt $\\nu$"          
                , 9:  "Cosmic"          
                , -1: "Unknown"
                , -99: "Bad Reco Signal"
            }

event_type = [ 0
              , 1
              #, 2
              , 3
              , 4
              , 5
              , 6
              , 7
              , 8
              , 9
              #, -99
              #, -1
              ]

# Open the file and get the tree
def file_globals(rootfile,file_type=""):
    file = uproot.open(rootfile)["hnlpizeroana"]
    events  = file["events"]
    subruns = file["subruns"]

    POTs = np.sum(subruns["pot"].array())
    if file_type=="cosmics": POTs = 1
    POTs_scaling =SBND_POTs/ POTs 
    N_SLICES = len(ak.flatten(events["slc_n_trks"].array()))
    
    if not file_type=="cosmics":sample_spills = np.sum(subruns['ngenevts'])
    else:                       sample_spills = 1
    
    spills = sample_spills * POTs_scaling

    #prepare labels for each event type
    event_types={}
    types=np.array(ak.flatten(events["slc_true_event_type"].array()))
    if file_type=="nu"     :    types[types==0 ] = -1
    if file_type=="cosmics":    types[types==-1] = 9
    for evt in event_type:
        event_types[evt] = types==evt

    globals_file={}
    globals_file["POTs_scaling"] = POTs_scaling
    globals_file["N_SLICES"] = N_SLICES
    globals_file["spills"] = spills
    globals_file["events"] = events
    globals_file["subruns"] = subruns
    globals_file["event_types"] = event_types
    
    # create a placeholder to store all the cuts, and update after each step
    globals_file["CUTS"] = np.ones(N_SLICES,dtype=bool)
    
    # in HNL case, first remove everything that does look like douplicated
    if file_type=="hnls": 
        ev_comp=ak.flatten(events['slc_comp']           .array() )>0.5
        nSig    =((types==0)*ev_comp)
        nNonFV  =((types==1)*ev_comp)

        globals_file["CUTS"] =np.array(nSig)
    
    return globals_file

def cosmics_reescale(g_cosmics,g_hnl,g_nus):

    target_pot = 1*10**21
    pot_per_spill = 5e12

    target_spill = target_pot / pot_per_spill 
    target_intime_spill = target_spill - g_hnl["spills"] - g_nus["spills"]
    sample_spill=np.sum(g_cosmics["subruns"]["ngenevts"].array())
    scale_pot = target_intime_spill / sample_spill
    g_cosmics["POTs_scaling"] = scale_pot


    print('-----------------------------------------------')
    print('target total spill = ' + str(target_spill))
    print('hnl + nu spill = ' + str(g_hnl["N_SLICES"] + g_nus["N_SLICES"]))
    print('target intime spill = ' + str(target_intime_spill))
    print('scale pot factor = ' + str(scale_pot))
    print('-----------------------------------------------')

    return g_cosmics

def print_efficiencies(g_hnl,g_nus,g_cosmics,nSig_t,nSig_r):
    
    hnl_sig  = sum(g_hnl["CUTS"]    *g_hnl["POTs_scaling"]*  (g_hnl["event_types"][0]))
    hnl_bckg = sum(g_hnl["CUTS"]    *g_hnl["POTs_scaling"]* ~(g_hnl["event_types"][0]))
    nus_bkg  = sum(g_nus["CUTS"]    *g_nus["POTs_scaling"])
    cos_bkg  = sum(g_cosmics["CUTS"]*g_cosmics["POTs_scaling"])

    # print("HNL signal: ",hnl_sig)
    # print("Nus       : ",nus_bkg)
    # print("Cosmics   : ",cos_bkg)
    eff = hnl_sig / nSig_t * 100
    select_eff = hnl_sig / nSig_r * 100
    purity = hnl_sig/(hnl_sig+nus_bkg+cos_bkg+hnl_bckg)*100

    print("---------------")
    print('purity = {0:.3g}'.format(purity))
    print('eff = {0:.3g}'.format(eff))
    print('select eff = {0:.3g}'.format(select_eff))

    return purity,eff,select_eff
