def _lepton_genpart_flav(event, lep_Hindex, lepton_position=2):
    """Return NanoAOD genPartFlav for one candidate-lepton position."""
    if lepton_position < 0 or lepton_position >= len(lep_Hindex):
        return 0
    lep_k = lep_Hindex[lepton_position]
    if lep_k < 0 or lep_k >= event.nLepton:
        return 0
    if lep_k < event.nElectron:
        return int(event.Electron_genPartFlav[lep_k])
    mu_k = lep_k - event.nElectron
    if mu_k < 0 or mu_k >= event.nMuon:
        return 0
    return int(event.Muon_genPartFlav[mu_k])


def PartOriginFromFlav(
    isData, event, lep_Hindex, FillValue, FillWeight, Hist_prompt,
    Hist_fakes, Hist_BDfakes, Hist_conv, isZX, lepton_position=2
):
    """
    Fill MC composition histograms from NanoAOD genPartFlav.

    ``lepton_position`` is the position inside the analysis candidate. 

    Mapping (Electron/Muon):
      1 -> prompt
      15 -> tau (counted in fakes,)
      22 -> conversion (electrons only; if present for muons, still treated as conv)
      5/4 -> heavy flavour (b/c) -> Hist_BDfakes
      3/0/other -> light or unmatched -> Hist_fakes
    """
    if isData:
        return

    flav = _lepton_genpart_flav(event, lep_Hindex, lepton_position)

    if flav == 1:
        if not isZX:
            Hist_prompt.Fill(FillValue, FillWeight)
    elif flav == 22:
        if not isZX:
            Hist_conv.Fill(FillValue, FillWeight)
    elif flav in (5, 4):
        if not isZX:
            Hist_BDfakes.Fill(FillValue, FillWeight)
    else:
        # Keep tau (15), light (3), unmatched (0), and unknown codes in fake bin
        if not isZX:
            Hist_fakes.Fill(FillValue, FillWeight)


def PartOrigin(isData, lep_matchedR03_PdgId, lep_matchedR03_MomId, lep_matchedR03_MomMomId, lep_Hindex, lep_id, FillValue, FillWeight, Hist_prompt, Hist_fakes, Hist_BDfakes,  Hist_conv, isZX):
    n_pass = 0
    n_nonprompttaus = 0
    n_nonpromptconv = 0
    n_nonpromptformphot = 0
    n_prompt = 0
    n_bdfake = 0
    n_fakes = 0
    n_others = 0  
    if not (isData):
        n_pass = n_pass + 1
        n=0
        test_particle_origin = ( (abs(lep_matchedR03_PdgId[lep_Hindex[2]])==15) or (abs(lep_matchedR03_MomId[lep_Hindex[2]])==15) or (abs(lep_matchedR03_MomId[lep_Hindex[2]])==15) )
        n = n + test_particle_origin
        
        
        if ( n >= 1 ):
            n_nonprompttaus = n_nonprompttaus + 1
            if not (isZX):
                Hist_fakes.Fill(FillValue,FillWeight)
            
        else:
            n = 0
            test_particle_origin = ( lep_matchedR03_PdgId[lep_Hindex[2]]==22 and (abs(lep_matchedR03_MomId[lep_Hindex[2]])==11 or abs(lep_matchedR03_MomId[lep_Hindex[2]])==13 or abs(lep_matchedR03_MomId[lep_Hindex[2]])==15) )
            n = n + test_particle_origin
            
            if (n==1): 
                n_nonpromptconv = n_nonpromptconv + 1
                if not (isZX):
                    Hist_conv.Fill(FillValue,FillWeight)

            else:
                n=0
                test_particle_origin = ( lep_matchedR03_PdgId[lep_Hindex[2]]==22 and not (abs(lep_matchedR03_MomId[lep_Hindex[2]])==11 or abs(lep_matchedR03_MomId[lep_Hindex[2]])==13 or abs(lep_matchedR03_MomId[lep_Hindex[2]])==15) )
                n = n + test_particle_origin
                
                # n+=(lep_matchedR03_PdgId[lep_Hindex[2]]==22)and((abs(lep_matchedR03_MomId[lep_Hindex[2]])==11)or(abs(lep_matchedR03_MomId[lep_Hindex[i]])==13));
                # n+=(lep_matchedR03_MomId[lep_Hindex[2]]==22)and((abs(lep_matchedR03_PdgId[lep_Hindex[2]])==11)or(abs(lep_matchedR03_PdgId[lep_Hindex[i]])==13));
                # n+=(lep_matchedR03_MomMomId[lep_Hindex[2]]==22)and((abs(lep_matchedR03_MomId[lep_Hindex[2]])==11)or(abs(lep_matchedR03_MomId[lep_Hindex[i]])==13));
                
                if (n >= 1):
                    if  not (isZX):
                        Hist_fakes.Fill(FillValue,FillWeight)
                    n_nonpromptformphot = n_nonpromptformphot + 1
        
                else:
                    n=0
                    
                    test_particle_origin_1 = ( lep_matchedR03_PdgId[lep_Hindex[2]]==23 )
                    test_particle_origin_2 = ( ( (lep_matchedR03_MomId[lep_Hindex[2]]==23) or (abs(lep_matchedR03_MomId[lep_Hindex[2]])==24) or (lep_matchedR03_MomId[lep_Hindex[2]]==25) ) and (lep_matchedR03_PdgId[lep_Hindex[2]]==lep_id[lep_Hindex[2]]) )
                    test_particle_origin_3 = ( ( (lep_matchedR03_MomMomId[lep_Hindex[2]]==23) or (abs(lep_matchedR03_MomMomId[lep_Hindex[2]])==24) or (lep_matchedR03_MomMomId[lep_Hindex[2]]==25) ) and ((lep_matchedR03_PdgId[lep_Hindex[2]]==lep_id[lep_Hindex[2]]) and (lep_matchedR03_MomId[lep_Hindex[2]]==lep_id[lep_Hindex[2]])))

                    n = n + test_particle_origin_1 + test_particle_origin_2 + test_particle_origin_3
                    if ( n >= 1 ):
                    
                        if not(isZX):
                            Hist_prompt.Fill(FillValue,FillWeight)
                        n_prompt = n_prompt + 1
                        
                    
                    else:
                        n=0

                        test_particle_origin_1 = ((abs(lep_matchedR03_PdgId[lep_Hindex[2]])>410)and(abs(lep_matchedR03_PdgId[lep_Hindex[2]])<436))or((abs(lep_matchedR03_PdgId[lep_Hindex[2]])>10410)and(abs(lep_matchedR03_PdgId[lep_Hindex[2]])<10433)) or ((abs(lep_matchedR03_PdgId[lep_Hindex[2]])>20412)and(abs(lep_matchedR03_PdgId[lep_Hindex[2]])<20434))
                        
                        test_particle_origin_2 = ((abs(lep_matchedR03_PdgId[lep_Hindex[2]])>510)and(abs(lep_matchedR03_PdgId[lep_Hindex[2]])<546))or((abs(lep_matchedR03_PdgId[lep_Hindex[2]])>10510)and(abs(lep_matchedR03_PdgId[lep_Hindex[2]])<10544)) or ((abs(lep_matchedR03_PdgId[lep_Hindex[2]])>20512)and(abs(lep_matchedR03_PdgId[lep_Hindex[2]])<20544))
                        
                        
                        test_particle_origin_3 = ((abs(lep_matchedR03_MomId[lep_Hindex[2]])>410)and(abs(lep_matchedR03_MomId[lep_Hindex[2]])<436))or((abs(lep_matchedR03_MomId[lep_Hindex[2]])>10410)and(abs(lep_matchedR03_MomId[lep_Hindex[2]])<10433)) or ((abs(lep_matchedR03_MomId[lep_Hindex[2]])>20412)and(abs(lep_matchedR03_MomId[lep_Hindex[2]])<20434))
                        
                        test_particle_origin_4 = ((abs(lep_matchedR03_MomId[lep_Hindex[2]])>510)and(abs(lep_matchedR03_MomId[lep_Hindex[2]])<546))or((abs(lep_matchedR03_MomId[lep_Hindex[2]])>10510)and(abs(abs(lep_matchedR03_MomId[lep_Hindex[2]]))<10544)) or ((abs(lep_matchedR03_MomId[lep_Hindex[2]])>20512)and(abs(lep_matchedR03_MomId[lep_Hindex[2]])<20544))
                        
                        
                        test_particle_origin_5 = ((abs(lep_matchedR03_MomMomId[lep_Hindex[2]])>410)and(abs(lep_matchedR03_MomMomId[lep_Hindex[2]])<436))or((abs(lep_matchedR03_MomMomId[lep_Hindex[2]])>10410)and(abs(lep_matchedR03_MomMomId[lep_Hindex[2]])<10433)) or ((abs(lep_matchedR03_MomMomId[lep_Hindex[2]])>20412)and(abs(lep_matchedR03_MomMomId[lep_Hindex[2]])<20434))
                        
                        test_particle_origin_6 = ((abs(lep_matchedR03_MomMomId[lep_Hindex[2]])>510)and(abs(lep_matchedR03_MomMomId[lep_Hindex[2]])<546))or((abs(lep_matchedR03_MomMomId[lep_Hindex[2]])>10510)and(abs(lep_matchedR03_MomMomId[lep_Hindex[2]])<10544)) or ((abs(lep_matchedR03_MomMomId[lep_Hindex[2]])>20512)and(abs(lep_matchedR03_MomMomId[lep_Hindex[2]])<20544))
                        
                        
                        n = n + test_particle_origin_1 + test_particle_origin_2 + test_particle_origin_3 + test_particle_origin_4 + test_particle_origin_5 + test_particle_origin_6
                        
                        if ( n >= 1 ) :
                            
                            if not(isZX):
                                Hist_BDfakes.Fill(FillValue,FillWeight)
                            
                            n_bdfake = n_bdfake + 1

                        else:
                            n = 0
                            
                            test_particle_origin_1 = ((abs(lep_matchedR03_PdgId[lep_Hindex[2]])<10)or(abs(lep_matchedR03_PdgId[lep_Hindex[2]])>100))
                            test_particle_origin_2 = ((abs(lep_matchedR03_MomId[lep_Hindex[2]])<10)or(abs(lep_matchedR03_MomId[lep_Hindex[2]])>100))
                            test_particle_origin_3 = ((abs(lep_matchedR03_MomMomId[lep_Hindex[2]])<10)or(abs(lep_matchedR03_MomMomId[lep_Hindex[2]])>100))

                            if ( n >= 1 ):
                                n_fakes = n_fakes + 1
                                
                                if not(isZX):
                                    Hist_fakes.Fill(FillValue,FillWeight)
                            
                            else:
                                n = 0
                                
                                if not(isZX):
                                    Hist_fakes.Fill(FillValue,FillWeight)
                                
                                n_others = n_others + 1
