import numpy as np
import awkward as ak

#templates for slc-lvl defined variables

def cut_var_equal(g,var,value):
    if type(g)==dict: # single file
        g["CUTS"]*= np.array(ak.flatten(g["events"][var].array())) == value
    if type(g)==tuple:# multiple files
        for gg in g:
            gg["CUTS"]*= np.array(ak.flatten(gg["events"][var].array())) ==value
    return g

def cut_var_more(g,var,value):
    if type(g)==dict: # single file
        g["CUTS"]*= np.array(ak.flatten(g["events"][var].array())) > value
    if type(g)==tuple:# multiple files
        for gg in g:
            gg["CUTS"]*= np.array(ak.flatten(gg["events"][var].array())) >value
    return g

def cut_var_less(g,var,value):
    if type(g)==dict: # single file
        g["CUTS"]*= np.array(ak.flatten(g["events"][var].array())) < value
    if type(g)==tuple:# multiple files
        for gg in g:
            gg["CUTS"]*= np.array(ak.flatten(gg["events"][var].array())) <value
    return g

def cut_pfp_track_like(g,score=0.6):
    if type(g)==tuple:# multiple files
        for gg in g:
            track_score=ak.flatten(gg["events"]["slc_pfp_track_score"].array())
            good=np.array(ak.sum(track_score<score,axis=1),dtype=bool)
            
            gg["CUTS"]*=good 
        return g
    else: raise Exception("Not implemented");

def cut_containement(g):
    #probably not the best way to do it
    if type(g)==tuple:# multiple files
        for gg in g:
            track      =ak.flatten(gg["events"]["slc_pfp_track_score"].array()) >=0.5
            track_end_x=ak.flatten(gg["events"]["slc_pfp_track_end_x"].array())
            track_end_y=ak.flatten(gg["events"]["slc_pfp_track_end_y"].array())
            track_end_z=ak.flatten(gg["events"]["slc_pfp_track_end_z"].array())
            
            shower=~track
            shower_end_x=ak.flatten(gg["events"]["slc_pfp_shower_end_x"].array())
            shower_end_y=ak.flatten(gg["events"]["slc_pfp_shower_end_y"].array())
            shower_end_z=ak.flatten(gg["events"]["slc_pfp_shower_end_z"].array())
            tracks_outside  = (abs(track_end_x )>195)|(abs(track_end_y  )>195) & (track_end_z<5) |(track_end_z>495)
            showers_outside = (abs(shower_end_x)>195)|(abs(shower_end_y )>195) & (shower_end_z<5)|(shower_end_z>495)
            tracks_outside=tracks_outside*track
            showers_outside=showers_outside*shower
            total=tracks_outside+showers_outside
            good=~np.array(ak.sum(total,axis=1),dtype=bool)
            gg["CUTS"]*=good 
        return total, good
        # return g
    else: raise Exception("Not implemented");

def cut_pfp_razzled_muon_score(g,score=0.08):
    if type(g)==tuple:# multiple files
        for gg in g:
            muon_score=ak.flatten(gg["events"]["slc_pfp_razzled_muon_score"].array())
            above=np.array(ak.sum(muon_score>score,axis=1),dtype=bool)
            
            gg["CUTS"]*=~above 
        return g
    else: raise Exception("Not implemented");

def cut_pfp_razzled_proton_score(g,score=0.4):
    if type(g)==tuple:# multiple files
        for gg in g:
            proton_score=ak.flatten(gg["events"]["slc_pfp_razzled_proton_score"].array())
            above=np.array(ak.sum(proton_score>score,axis=1),dtype=bool)
            
            gg["CUTS"]*=~above 
        return g
    else: raise Exception("Not implemented");

def cut_pfp_razzled_pion_score(g,score=0.8):
    if type(g)==tuple:# multiple files
        for gg in g:
            pion_score=ak.flatten(gg["events"]["slc_pfp_razzled_pion_score"].array())
            above=np.array(ak.sum(pion_score>score,axis=1),dtype=bool)
            
            gg["CUTS"]*=~above 
        return g
    else: raise Exception("Not implemented");
