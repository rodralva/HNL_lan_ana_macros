# %%
import uproot
import numpy as np
import awkward as ak
import matplotlib.pyplot as plt
from aux_functions import *
from cut_functions import *
import math

SBND_POTs   = 1e21
DEBUG=False

# %%
globals_hnl     = file_globals("../hnl_overlay_eenu_no_weightsosmics_2e-5.root",file_type="hnls")
globals_nus     = file_globals("../nu_5k.root"     ,file_type="nu")
globals_cosmics = file_globals("../cosmics_5k.root",file_type="cosmics")
cosmics_reescale(globals_cosmics,globals_hnl,globals_nus)
if DEBUG:
    print("POTs scaling HNL: ",globals_hnl["POTs_scaling"])
    print("POTs scaling NUS: ",globals_nus["POTs_scaling"])
    print("POTs scaling COS: ",globals_cosmics["POTs_scaling"])

nSig_t, nNonFV_t=get_true_signal_in_all_spills(globals_hnl["events"],globals_hnl["POTs_scaling"],DEBUG=True)
nSig_r, nNonFV_r=get_reco_signal_in_all_spills(globals_hnl["events"],globals_hnl["POTs_scaling"],DEBUG=True)

print_efficiencies(globals_hnl,globals_nus,globals_cosmics,nSig_t,nSig_r);
print("CUTS: ","noCut")

#CUT presel
cut_var_equal((globals_hnl,globals_nus,globals_cosmics),"slc_is_clear_cosmics",0)
cut_var_equal((globals_hnl,globals_nus,globals_cosmics),"slc_is_fv"           ,1)
cut_pfp_track_like((globals_hnl,globals_nus,globals_cosmics))
print_efficiencies(globals_hnl,globals_nus,globals_cosmics,nSig_t,nSig_r);
print("CUTS: ","presel")

#COSMICS
cut_var_more((globals_hnl,globals_nus,globals_cosmics),"slc_crumbs_score"      ,0)
print_efficiencies(globals_hnl,globals_nus,globals_cosmics,nSig_t,nSig_r);
print("CUTS: ","cosmic")

#CUT_between_buckets
cutBetweenBucket((globals_nus,globals_cosmics,globals_hnl));
print_efficiencies(globals_hnl,globals_nus,globals_cosmics,nSig_t,nSig_r);
print("CUTS: ","BetweenBuckets")

#CONTAINEMENT
total, good=cut_containement((globals_hnl,globals_nus,globals_cosmics))
print_efficiencies(globals_hnl,globals_nus,globals_cosmics,nSig_t,nSig_r);
print("CUTS: ","contained")

#muon PID
cut_var_equal((globals_hnl,globals_nus,globals_cosmics),"slc_n_razzled_muons",0);
cut_pfp_razzled_muon_score((globals_hnl,globals_nus,globals_cosmics));
print_efficiencies(globals_hnl,globals_nus,globals_cosmics,nSig_t,nSig_r);
print("CUTS: ","muon")

#proton PID
cut_var_less((globals_hnl,globals_nus,globals_cosmics),"slc_n_razzled_protons_thresh",2);
cut_pfp_razzled_proton_score((globals_hnl,globals_nus,globals_cosmics));
print_efficiencies(globals_hnl,globals_nus,globals_cosmics,nSig_t,nSig_r);
print("CUTS: ","proton")

#pion PID
cut_var_less((globals_hnl,globals_nus,globals_cosmics),"slc_n_razzled_pions_thresh",2);
cut_pfp_razzled_pion_score((globals_hnl,globals_nus,globals_cosmics));
print_efficiencies(globals_hnl,globals_nus,globals_cosmics,nSig_t,nSig_r);
print("CUTS: ","pion")

#shower req
cut_var_more((globals_hnl,globals_nus,globals_cosmics),"slc_n_shws",0);
cut_var_less((globals_hnl,globals_nus,globals_cosmics),"slc_n_shws",5);
cut_shower_angle((globals_hnl,globals_nus,globals_cosmics),value=30)
print_efficiencies(globals_hnl,globals_nus,globals_cosmics,nSig_t,nSig_r);
print("CUTS: ","shower")



