//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Sep 16 16:24:03 2016 by ROOT version 5.34/36
// from TTree tree/H2TauTauSyncTreeProducerTauMu
// found on file: /afs/cern.ch/user/s/steggema/public/Sync2016/mt_v2.root
//////////////////////////////////////////////////////////

#ifndef SynchCERN_h
#define SynchCERN_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class SynchCERN {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain
   std::string 	   suffix_;  //name of the sample (e.g. WAW, CERN, KIT, Vienna)
   
   //pileup
   TH1F *hNPV;
   TH1F *hNPU;
   //muon
   TH1F*hMuPt;
   TH1F*hMuPhi;
   TH1F*hMuEta;
   TH1F*hMud0;
   TH1F*hMudZ;
   TH1F*hMuIsoZoom;
   TH1F*hMuIso;
   TH1F*hMuMt;
   TH1F*hMuPFMt;
   TH1F*hMuM;
   TH1F*hMuGenMatch;
   //tau
   TH1F*hTauPt;
   TH1F*hTauPhi;
   TH1F*hTauEta;
   TH1F*hTaud0;
   TH1F*hTaudZ;
   TH1F*hTauIso;
   TH1F*hTauMt;
   TH1F*hTauPFMt;
   TH1F*hTauGenMatch;
   //Jets
   TH1F*hJet1Pt;
   TH1F*hJet1Phi;
   TH1F*hJet1Eta;
   TH1F*hJet2Pt;
   TH1F*hJet2Phi;
   TH1F*hJet2Eta;
   TH1F*hBJet1Pt;
   TH1F*hBJet1Phi;
   TH1F*hBJet1Eta;
   TH1F*hBJet2Pt;
   TH1F*hBJet2Phi;
   TH1F*hBJet2Eta;
   //MET
   TH1F*hMET;
   TH1F*hMETphi;
   TH1F*hMVAMET;
   TH1F*hMVAMETphi;
   TH1F*hMVACov00;
   TH1F*hMVACov01;
   TH1F*hMVACov10;
   TH1F*hMVACov11;
   TH1F*hPZetaVis;
   TH1F*hPZetaMiss;
   TH1F*hPFPZetaMiss;
   //di-tau
   TH1F*hPtTT;
   TH1F*hMtTot;
   TH1F*hMVis;
   TH1F*hMSV;
   //vbf system
   TH1F *hMJJ;
   TH1F *hJdEta;
   TH1F *hJdPhi;
   //extra lepton vetos
   TH1F *hDiLeptonVeto;
   TH1F *hExtraMuonVeto;
   TH1F *hExtraElectronVeto;

   // Declaration of leaf types
   Int_t           run;
   Int_t           lumi;
   Int_t           evt;
   Int_t           bx;
   Int_t           orbit_number;
   Int_t           is_data;
   Double_t        npu;
   Int_t           pass_leptons;
   Int_t           dilepton_veto;
   Int_t           extramuon_veto;
   Int_t           extraelec_veto;
   Int_t           njets;
   Int_t           n_jets_puid;
   Int_t           njetspt20;
   Int_t           n_jets_20_puid;
   Int_t           nbtag;
   Int_t           n_jets_csvl;
   Int_t           npv;
   Double_t        rho;
   Double_t        weight;
   Double_t        puweight;
   Double_t        weight_njet;
   Double_t        m_vis;
   Double_t        mt_total;
   Double_t        sum_lepton_mt;
   Double_t        sqsum_lepton_mt;
   Double_t        pzetamiss;
   Double_t        pzetavis;
   Double_t        pzeta_disc;
   Double_t        mt_1;
   Double_t        mt_2;
   Double_t        mvacov00;
   Double_t        mvacov10;
   Double_t        mvacov11;
   Double_t        mvametphi;
   Double_t        met_px;
   Double_t        met_py;
   Double_t        mvamet;
   Double_t        pt_tt;
   Double_t        delta_phi_l1_l2;
   Double_t        delta_eta_l1_l2;
   Double_t        delta_r_l1_l2;
   Double_t        delta_phi_l1_met;
   Double_t        delta_phi_l2_met;
   Double_t        m_sv;
   Double_t        mt_sv;
   Double_t        svfit_mass_error;
   Double_t        pt_sv;
   Double_t        svfit_l1_pt;
   Double_t        svfit_l1_eta;
   Double_t        svfit_l1_phi;
   Double_t        svfit_l1_charge;
   Double_t        svfit_l1_mass;
   Double_t        svfit_l2_pt;
   Double_t        svfit_l2_eta;
   Double_t        svfit_l2_phi;
   Double_t        svfit_l2_charge;
   Double_t        svfit_l2_mass;
   Double_t        geninfo_mcweight;
   Int_t           NUP;
   Double_t        geninfo_htgen;
   Double_t        geninvmass;
   Double_t        weight_gen;
   Double_t        genmet_pt;
   Double_t        genmet_px;
   Double_t        genmet_py;
   Double_t        genmet_phi;
   Double_t        vbf_mjj;
   Double_t        vbf_deta;
   Int_t           vbf_n_central20;
   Int_t           vbf_n_central;
   Double_t        vbf_jdphi;
   Double_t        vbf_dijetpt;
   Double_t        vbf_dijetphi;
   Double_t        vbf_dphidijethiggs;
   Double_t        vbf_mindetajetvis;
   Double_t        jpt_1;
   Double_t        jeta_1;
   Double_t        jphi_1;
   Double_t        jet1_charge;
   Double_t        jet1_mass;
   Double_t        jmva_1;
   Double_t        jpuid_1;
   Double_t        jet1_flavour_parton;
   Double_t        jcsv_1;
   Double_t        jet1_genjet_pt;
   Double_t        jpt_2;
   Double_t        jeta_2;
   Double_t        jphi_2;
   Double_t        jet2_charge;
   Double_t        jet2_mass;
   Double_t        jmva_2;
   Double_t        jpuid_2;
   Double_t        jet2_flavour_parton;
   Double_t        jcsv_2;
   Double_t        jet2_genjet_pt;
   Double_t        bpt_1;
   Double_t        beta_1;
   Double_t        bphi_1;
   Double_t        bjet1_charge;
   Double_t        bjet1_mass;
   Double_t        bmva_1;
   Double_t        bpuid_1;
   Double_t        bjet1_flavour_parton;
   Double_t        bcsv_1;
   Double_t        bjet1_genjet_pt;
   Double_t        bpt_2;
   Double_t        beta_2;
   Double_t        bphi_2;
   Double_t        bjet2_charge;
   Double_t        bjet2_mass;
   Double_t        bmva_2;
   Double_t        bpuid_2;
   Double_t        bjet2_flavour_parton;
   Double_t        bcsv_2;
   Double_t        bjet2_genjet_pt;
   Double_t        HT_allJets;
   Double_t        HT_jets;
   Double_t        HT_bJets;
   Double_t        HT_cleanJets;
   Double_t        HT_jets30;
   Double_t        HT_cleanJets30;
   Double_t        genboson_pt;
   Double_t        genboson_eta;
   Double_t        genboson_phi;
   Double_t        genboson_charge;
   Double_t        genboson_mass;
   Double_t        genboson_pdgId;
   Double_t        gen_top_1_pt;
   Double_t        gen_top_2_pt;
   Double_t        gen_top_weight;
   Double_t        puppimet;
   Double_t        puppimetphi;
   Double_t        puppimt_1;
   Double_t        puppimt_2;
   Double_t        met;
   Double_t        metphi;
   Double_t        pfmet_mt1;
   Double_t        pfmet_mt2;
   Double_t        pt_2;
   Double_t        eta_2;
   Double_t        phi_2;
   Double_t        q_2;
   Double_t        m_2;
   Double_t        l2_jet_pt;
   Double_t        l2_jet_eta;
   Double_t        l2_jet_phi;
   Double_t        l2_jet_charge;
   Double_t        l2_jet_mass;
   Double_t        d0_2;
   Double_t        l2_dxy_error;
   Double_t        dZ_2;
   Double_t        l2_dz_error;
   Double_t        l2_weight;
   Double_t        trigweight_2;
   Double_t        l2_weight_eff_data_trigger;
   Double_t        l2_eff_trigger_data;
   Double_t        l2_eff_trigger_mc;
   Double_t        l2_weight_idiso;
   Double_t        l2_eff_idiso_data;
   Double_t        l2_eff_idiso_mc;
   Double_t        gen_match_2;
   Double_t        l2_decayMode;
   Double_t        l2_zImpact;
   Double_t        l2_dz_selfvertex;
   Double_t        l2_ptScale;
   Int_t           l2_againstElectronMVA6;
   Int_t           l2_againstMuon3;
   Double_t        byCombinedIsolationDeltaBetaCorrRaw3Hits_2;
   Double_t        l2_byIsolationMVArun2v1DBoldDMwLTraw;
   Double_t        l2_byIsolationMVArun2v1DBnewDMwLTraw;
   Double_t        l2_byIsolationMVArun2v1DBdR03oldDMwLTraw;
   Int_t           l2_byCombinedIsolationDeltaBetaCorr3Hits;
   Int_t           l2_byIsolationMVArun2v1DBoldDMwLT;
   Int_t           l2_byIsolationMVArun2v1DBnewDMwLT;
   Int_t           l2_byIsolationMVArun2v1DBdR03oldDMwLT;
   Double_t        chargedIsoPtSum_2;
   Double_t        decayModeFindingOldDMs_2;
   Double_t        l2_footprintCorrection;
   Double_t        neutralIsoPtSum_2;
   Double_t        puCorrPtSum_2;
   Double_t        l2_photonPtSumOutsideSignalCone;
   Double_t        mva_olddm_tight_2;
   Double_t        pt_1;
   Double_t        eta_1;
   Double_t        phi_1;
   Double_t        q_1;
   Double_t        m_1;
   Double_t        l1_jet_pt;
   Double_t        l1_jet_eta;
   Double_t        l1_jet_phi;
   Double_t        l1_jet_charge;
   Double_t        l1_jet_mass;
   Double_t        d0_1;
   Double_t        l1_dxy_error;
   Double_t        dZ_1;
   Double_t        l1_dz_error;
   Double_t        l1_weight;
   Double_t        trigweight_1;
   Double_t        l1_weight_eff_data_trigger;
   Double_t        l1_eff_trigger_data;
   Double_t        l1_eff_trigger_mc;
   Double_t        idisoweight_1;
   Double_t        l1_eff_idiso_data;
   Double_t        l1_eff_idiso_mc;
   Double_t        gen_match_1;
   Double_t        iso_1;
   Double_t        l1_reliso05_04;
   Double_t        id_m_loose_1;
   Double_t        id_m_medium_1;
   Double_t        id_m_tight_1;
   Double_t        id_m_tightnovtx_1;
   Double_t        id_m_highpt_1;
   Double_t        l1_dxy_innertrack;
   Double_t        l1_dz_innertrack;
   Double_t        l2_gen_pt;
   Double_t        l2_gen_eta;
   Double_t        l2_gen_phi;
   Double_t        l2_gen_charge;
   Double_t        l2_gen_mass;
   Double_t        l2_gen_pdgId;
   Int_t           l2_gen_lepfromtau;
   Double_t        l1_gen_pt;
   Double_t        l1_gen_eta;
   Double_t        l1_gen_phi;
   Double_t        l1_gen_charge;
   Double_t        l1_gen_mass;
   Double_t        l1_gen_pdgId;
   Int_t           l1_gen_lepfromtau;
   Double_t        l2_gen_vis_pt;
   Double_t        l2_gen_vis_eta;
   Double_t        l2_gen_vis_phi;
   Double_t        l2_gen_vis_charge;
   Double_t        l2_gen_vis_mass;
   Int_t           l2_gen_decaymode;
   Double_t        l2_gen_nc_ratio;
   Double_t        l2_nc_ratio;
   Double_t        l2_weight_fakerate;
   Double_t        l2_weight_fakerate_up;
   Double_t        l2_weight_fakerate_down;

   // List of branches
   TBranch        *b_run;   //!
   TBranch        *b_lumi;   //!
   TBranch        *b_evt;   //!
   TBranch        *b_bx;   //!
   TBranch        *b_orbit_number;   //!
   TBranch        *b_is_data;   //!
   TBranch        *b_npu;   //!
   TBranch        *b_pass_leptons;   //!
   TBranch        *b_dilepton_veto;   //!
   TBranch        *b_extramuon_veto;   //!
   TBranch        *b_extraelec_veto;   //!
   TBranch        *b_njets;   //!
   TBranch        *b_n_jets_puid;   //!
   TBranch        *b_njetspt20;   //!
   TBranch        *b_n_jets_20_puid;   //!
   TBranch        *b_nbtag;   //!
   TBranch        *b_n_jets_csvl;   //!
   TBranch        *b_npv;   //!
   TBranch        *b_rho;   //!
   TBranch        *b_weight;   //!
   TBranch        *b_puweight;   //!
   TBranch        *b_weight_njet;   //!
   TBranch        *b_m_vis;   //!
   TBranch        *b_mt_total;   //!
   TBranch        *b_sum_lepton_mt;   //!
   TBranch        *b_sqsum_lepton_mt;   //!
   TBranch        *b_pzetamiss;   //!
   TBranch        *b_pzetavis;   //!
   TBranch        *b_pzeta_disc;   //!
   TBranch        *b_mt_1;   //!
   TBranch        *b_mt_2;   //!
   TBranch        *b_mvacov00;   //!
   TBranch        *b_mvacov10;   //!
   TBranch        *b_mvacov11;   //!
   TBranch        *b_mvametphi;   //!
   TBranch        *b_met_px;   //!
   TBranch        *b_met_py;   //!
   TBranch        *b_mvamet;   //!
   TBranch        *b_pt_tt;   //!
   TBranch        *b_delta_phi_l1_l2;   //!
   TBranch        *b_delta_eta_l1_l2;   //!
   TBranch        *b_delta_r_l1_l2;   //!
   TBranch        *b_delta_phi_l1_met;   //!
   TBranch        *b_delta_phi_l2_met;   //!
   TBranch        *b_m_sv;   //!
   TBranch        *b_mt_sv;   //!
   TBranch        *b_svfit_mass_error;   //!
   TBranch        *b_pt_sv;   //!
   TBranch        *b_svfit_l1_pt;   //!
   TBranch        *b_svfit_l1_eta;   //!
   TBranch        *b_svfit_l1_phi;   //!
   TBranch        *b_svfit_l1_charge;   //!
   TBranch        *b_svfit_l1_mass;   //!
   TBranch        *b_svfit_l2_pt;   //!
   TBranch        *b_svfit_l2_eta;   //!
   TBranch        *b_svfit_l2_phi;   //!
   TBranch        *b_svfit_l2_charge;   //!
   TBranch        *b_svfit_l2_mass;   //!
   TBranch        *b_geninfo_mcweight;   //!
   TBranch        *b_NUP;   //!
   TBranch        *b_geninfo_htgen;   //!
   TBranch        *b_geninvmass;   //!
   TBranch        *b_weight_gen;   //!
   TBranch        *b_genmet_pt;   //!
   TBranch        *b_genmet_px;   //!
   TBranch        *b_genmet_py;   //!
   TBranch        *b_genmet_phi;   //!
   TBranch        *b_vbf_mjj;   //!
   TBranch        *b_vbf_deta;   //!
   TBranch        *b_vbf_n_central20;   //!
   TBranch        *b_vbf_n_central;   //!
   TBranch        *b_vbf_jdphi;   //!
   TBranch        *b_vbf_dijetpt;   //!
   TBranch        *b_vbf_dijetphi;   //!
   TBranch        *b_vbf_dphidijethiggs;   //!
   TBranch        *b_vbf_mindetajetvis;   //!
   TBranch        *b_jpt_1;   //!
   TBranch        *b_jeta_1;   //!
   TBranch        *b_jphi_1;   //!
   TBranch        *b_jet1_charge;   //!
   TBranch        *b_jet1_mass;   //!
   TBranch        *b_jmva_1;   //!
   TBranch        *b_jpuid_1;   //!
   TBranch        *b_jet1_flavour_parton;   //!
   TBranch        *b_jcsv_1;   //!
   TBranch        *b_jet1_genjet_pt;   //!
   TBranch        *b_jpt_2;   //!
   TBranch        *b_jeta_2;   //!
   TBranch        *b_jphi_2;   //!
   TBranch        *b_jet2_charge;   //!
   TBranch        *b_jet2_mass;   //!
   TBranch        *b_jmva_2;   //!
   TBranch        *b_jpuid_2;   //!
   TBranch        *b_jet2_flavour_parton;   //!
   TBranch        *b_jcsv_2;   //!
   TBranch        *b_jet2_genjet_pt;   //!
   TBranch        *b_bpt_1;   //!
   TBranch        *b_beta_1;   //!
   TBranch        *b_bphi_1;   //!
   TBranch        *b_bjet1_charge;   //!
   TBranch        *b_bjet1_mass;   //!
   TBranch        *b_bmva_1;   //!
   TBranch        *b_bpuid_1;   //!
   TBranch        *b_bjet1_flavour_parton;   //!
   TBranch        *b_bcsv_1;   //!
   TBranch        *b_bjet1_genjet_pt;   //!
   TBranch        *b_bpt_2;   //!
   TBranch        *b_beta_2;   //!
   TBranch        *b_bphi_2;   //!
   TBranch        *b_bjet2_charge;   //!
   TBranch        *b_bjet2_mass;   //!
   TBranch        *b_bmva_2;   //!
   TBranch        *b_bpuid_2;   //!
   TBranch        *b_bjet2_flavour_parton;   //!
   TBranch        *b_bcsv_2;   //!
   TBranch        *b_bjet2_genjet_pt;   //!
   TBranch        *b_HT_allJets;   //!
   TBranch        *b_HT_jets;   //!
   TBranch        *b_HT_bJets;   //!
   TBranch        *b_HT_cleanJets;   //!
   TBranch        *b_HT_jets30;   //!
   TBranch        *b_HT_cleanJets30;   //!
   TBranch        *b_genboson_pt;   //!
   TBranch        *b_genboson_eta;   //!
   TBranch        *b_genboson_phi;   //!
   TBranch        *b_genboson_charge;   //!
   TBranch        *b_genboson_mass;   //!
   TBranch        *b_genboson_pdgId;   //!
   TBranch        *b_gen_top_1_pt;   //!
   TBranch        *b_gen_top_2_pt;   //!
   TBranch        *b_gen_top_weight;   //!
   TBranch        *b_puppimet;   //!
   TBranch        *b_puppimetphi;   //!
   TBranch        *b_puppimt_1;   //!
   TBranch        *b_puppimt_2;   //!
   TBranch        *b_met;   //!
   TBranch        *b_metphi;   //!
   TBranch        *b_pfmet_mt1;   //!
   TBranch        *b_pfmet_mt2;   //!
   TBranch        *b_pt_2;   //!
   TBranch        *b_eta_2;   //!
   TBranch        *b_phi_2;   //!
   TBranch        *b_q_2;   //!
   TBranch        *b_m_2;   //!
   TBranch        *b_l2_jet_pt;   //!
   TBranch        *b_l2_jet_eta;   //!
   TBranch        *b_l2_jet_phi;   //!
   TBranch        *b_l2_jet_charge;   //!
   TBranch        *b_l2_jet_mass;   //!
   TBranch        *b_d0_2;   //!
   TBranch        *b_l2_dxy_error;   //!
   TBranch        *b_dZ_2;   //!
   TBranch        *b_l2_dz_error;   //!
   TBranch        *b_l2_weight;   //!
   TBranch        *b_trigweight_2;   //!
   TBranch        *b_l2_weight_eff_data_trigger;   //!
   TBranch        *b_l2_eff_trigger_data;   //!
   TBranch        *b_l2_eff_trigger_mc;   //!
   TBranch        *b_l2_weight_idiso;   //!
   TBranch        *b_l2_eff_idiso_data;   //!
   TBranch        *b_l2_eff_idiso_mc;   //!
   TBranch        *b_gen_match_2;   //!
   TBranch        *b_l2_decayMode;   //!
   TBranch        *b_l2_zImpact;   //!
   TBranch        *b_l2_dz_selfvertex;   //!
   TBranch        *b_l2_ptScale;   //!
   TBranch        *b_l2_againstElectronMVA6;   //!
   TBranch        *b_l2_againstMuon3;   //!
   TBranch        *b_byCombinedIsolationDeltaBetaCorrRaw3Hits_2;   //!
   TBranch        *b_l2_byIsolationMVArun2v1DBoldDMwLTraw;   //!
   TBranch        *b_l2_byIsolationMVArun2v1DBnewDMwLTraw;   //!
   TBranch        *b_l2_byIsolationMVArun2v1DBdR03oldDMwLTraw;   //!
   TBranch        *b_l2_byCombinedIsolationDeltaBetaCorr3Hits;   //!
   TBranch        *b_l2_byIsolationMVArun2v1DBoldDMwLT;   //!
   TBranch        *b_l2_byIsolationMVArun2v1DBnewDMwLT;   //!
   TBranch        *b_l2_byIsolationMVArun2v1DBdR03oldDMwLT;   //!
   TBranch        *b_chargedIsoPtSum_2;   //!
   TBranch        *b_decayModeFindingOldDMs_2;   //!
   TBranch        *b_l2_footprintCorrection;   //!
   TBranch        *b_neutralIsoPtSum_2;   //!
   TBranch        *b_puCorrPtSum_2;   //!
   TBranch        *b_l2_photonPtSumOutsideSignalCone;   //!
   TBranch        *b_mva_olddm_tight_2;   //!
   TBranch        *b_pt_1;   //!
   TBranch        *b_eta_1;   //!
   TBranch        *b_phi_1;   //!
   TBranch        *b_q_1;   //!
   TBranch        *b_m_1;   //!
   TBranch        *b_l1_jet_pt;   //!
   TBranch        *b_l1_jet_eta;   //!
   TBranch        *b_l1_jet_phi;   //!
   TBranch        *b_l1_jet_charge;   //!
   TBranch        *b_l1_jet_mass;   //!
   TBranch        *b_d0_1;   //!
   TBranch        *b_l1_dxy_error;   //!
   TBranch        *b_dZ_1;   //!
   TBranch        *b_l1_dz_error;   //!
   TBranch        *b_l1_weight;   //!
   TBranch        *b_trigweight_1;   //!
   TBranch        *b_l1_weight_eff_data_trigger;   //!
   TBranch        *b_l1_eff_trigger_data;   //!
   TBranch        *b_l1_eff_trigger_mc;   //!
   TBranch        *b_idisoweight_1;   //!
   TBranch        *b_l1_eff_idiso_data;   //!
   TBranch        *b_l1_eff_idiso_mc;   //!
   TBranch        *b_gen_match_1;   //!
   TBranch        *b_iso_1;   //!
   TBranch        *b_l1_reliso05_04;   //!
   TBranch        *b_id_m_loose_1;   //!
   TBranch        *b_id_m_medium_1;   //!
   TBranch        *b_id_m_tight_1;   //!
   TBranch        *b_id_m_tightnovtx_1;   //!
   TBranch        *b_id_m_highpt_1;   //!
   TBranch        *b_l1_dxy_innertrack;   //!
   TBranch        *b_l1_dz_innertrack;   //!
   TBranch        *b_l2_gen_pt;   //!
   TBranch        *b_l2_gen_eta;   //!
   TBranch        *b_l2_gen_phi;   //!
   TBranch        *b_l2_gen_charge;   //!
   TBranch        *b_l2_gen_mass;   //!
   TBranch        *b_l2_gen_pdgId;   //!
   TBranch        *b_l2_gen_lepfromtau;   //!
   TBranch        *b_l1_gen_pt;   //!
   TBranch        *b_l1_gen_eta;   //!
   TBranch        *b_l1_gen_phi;   //!
   TBranch        *b_l1_gen_charge;   //!
   TBranch        *b_l1_gen_mass;   //!
   TBranch        *b_l1_gen_pdgId;   //!
   TBranch        *b_l1_gen_lepfromtau;   //!
   TBranch        *b_l2_gen_vis_pt;   //!
   TBranch        *b_l2_gen_vis_eta;   //!
   TBranch        *b_l2_gen_vis_phi;   //!
   TBranch        *b_l2_gen_vis_charge;   //!
   TBranch        *b_l2_gen_vis_mass;   //!
   TBranch        *b_l2_gen_decaymode;   //!
   TBranch        *b_l2_gen_nc_ratio;   //!
   TBranch        *b_l2_nc_ratio;   //!
   TBranch        *b_l2_weight_fakerate;   //!
   TBranch        *b_l2_weight_fakerate_up;   //!
   TBranch        *b_l2_weight_fakerate_down;   //!

   SynchCERN(TTree *tree=0, std::string suffix = "");
   virtual ~SynchCERN();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef SynchCERN_cxx
SynchCERN::SynchCERN(TTree *tree, std::string suffix) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   suffix_ = suffix;
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/afs/cern.ch/user/s/steggema/public/Sync2016/mt_v2.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("/afs/cern.ch/user/s/steggema/public/Sync2016/mt_v2.root");
      }
      f->GetObject("tree",tree);

   }
   Init(tree);
}

SynchCERN::~SynchCERN()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t SynchCERN::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t SynchCERN::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void SynchCERN::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("lumi", &lumi, &b_lumi);
   fChain->SetBranchAddress("evt", &evt, &b_evt);
   fChain->SetBranchAddress("bx", &bx, &b_bx);
   fChain->SetBranchAddress("orbit_number", &orbit_number, &b_orbit_number);
   fChain->SetBranchAddress("is_data", &is_data, &b_is_data);
   fChain->SetBranchAddress("npu", &npu, &b_npu);
   fChain->SetBranchAddress("pass_leptons", &pass_leptons, &b_pass_leptons);
   fChain->SetBranchAddress("dilepton_veto", &dilepton_veto, &b_dilepton_veto);
   fChain->SetBranchAddress("extramuon_veto", &extramuon_veto, &b_extramuon_veto);
   fChain->SetBranchAddress("extraelec_veto", &extraelec_veto, &b_extraelec_veto);
   fChain->SetBranchAddress("njets", &njets, &b_njets);
   fChain->SetBranchAddress("n_jets_puid", &n_jets_puid, &b_n_jets_puid);
   fChain->SetBranchAddress("njetspt20", &njetspt20, &b_njetspt20);
   fChain->SetBranchAddress("n_jets_20_puid", &n_jets_20_puid, &b_n_jets_20_puid);
   fChain->SetBranchAddress("nbtag", &nbtag, &b_nbtag);
   fChain->SetBranchAddress("n_jets_csvl", &n_jets_csvl, &b_n_jets_csvl);
   fChain->SetBranchAddress("npv", &npv, &b_npv);
   fChain->SetBranchAddress("rho", &rho, &b_rho);
   fChain->SetBranchAddress("weight", &weight, &b_weight);
   fChain->SetBranchAddress("puweight", &puweight, &b_puweight);
   fChain->SetBranchAddress("weight_njet", &weight_njet, &b_weight_njet);
   fChain->SetBranchAddress("m_vis", &m_vis, &b_m_vis);
   fChain->SetBranchAddress("mt_total", &mt_total, &b_mt_total);
   fChain->SetBranchAddress("sum_lepton_mt", &sum_lepton_mt, &b_sum_lepton_mt);
   fChain->SetBranchAddress("sqsum_lepton_mt", &sqsum_lepton_mt, &b_sqsum_lepton_mt);
   fChain->SetBranchAddress("pzetamiss", &pzetamiss, &b_pzetamiss);
   fChain->SetBranchAddress("pzetavis", &pzetavis, &b_pzetavis);
   fChain->SetBranchAddress("pzeta_disc", &pzeta_disc, &b_pzeta_disc);
   fChain->SetBranchAddress("mt_1", &mt_1, &b_mt_1);
   fChain->SetBranchAddress("mt_2", &mt_2, &b_mt_2);
   fChain->SetBranchAddress("mvacov00", &mvacov00, &b_mvacov00);
   fChain->SetBranchAddress("mvacov10", &mvacov10, &b_mvacov10);
   fChain->SetBranchAddress("mvacov11", &mvacov11, &b_mvacov11);
   fChain->SetBranchAddress("mvametphi", &mvametphi, &b_mvametphi);
   fChain->SetBranchAddress("met_px", &met_px, &b_met_px);
   fChain->SetBranchAddress("met_py", &met_py, &b_met_py);
   fChain->SetBranchAddress("mvamet", &mvamet, &b_mvamet);
   fChain->SetBranchAddress("pt_tt", &pt_tt, &b_pt_tt);
   fChain->SetBranchAddress("delta_phi_l1_l2", &delta_phi_l1_l2, &b_delta_phi_l1_l2);
   fChain->SetBranchAddress("delta_eta_l1_l2", &delta_eta_l1_l2, &b_delta_eta_l1_l2);
   fChain->SetBranchAddress("delta_r_l1_l2", &delta_r_l1_l2, &b_delta_r_l1_l2);
   fChain->SetBranchAddress("delta_phi_l1_met", &delta_phi_l1_met, &b_delta_phi_l1_met);
   fChain->SetBranchAddress("delta_phi_l2_met", &delta_phi_l2_met, &b_delta_phi_l2_met);
   fChain->SetBranchAddress("m_sv", &m_sv, &b_m_sv);
   fChain->SetBranchAddress("mt_sv", &mt_sv, &b_mt_sv);
   fChain->SetBranchAddress("svfit_mass_error", &svfit_mass_error, &b_svfit_mass_error);
   fChain->SetBranchAddress("pt_sv", &pt_sv, &b_pt_sv);
   fChain->SetBranchAddress("svfit_l1_pt", &svfit_l1_pt, &b_svfit_l1_pt);
   fChain->SetBranchAddress("svfit_l1_eta", &svfit_l1_eta, &b_svfit_l1_eta);
   fChain->SetBranchAddress("svfit_l1_phi", &svfit_l1_phi, &b_svfit_l1_phi);
   fChain->SetBranchAddress("svfit_l1_charge", &svfit_l1_charge, &b_svfit_l1_charge);
   fChain->SetBranchAddress("svfit_l1_mass", &svfit_l1_mass, &b_svfit_l1_mass);
   fChain->SetBranchAddress("svfit_l2_pt", &svfit_l2_pt, &b_svfit_l2_pt);
   fChain->SetBranchAddress("svfit_l2_eta", &svfit_l2_eta, &b_svfit_l2_eta);
   fChain->SetBranchAddress("svfit_l2_phi", &svfit_l2_phi, &b_svfit_l2_phi);
   fChain->SetBranchAddress("svfit_l2_charge", &svfit_l2_charge, &b_svfit_l2_charge);
   fChain->SetBranchAddress("svfit_l2_mass", &svfit_l2_mass, &b_svfit_l2_mass);
   fChain->SetBranchAddress("geninfo_mcweight", &geninfo_mcweight, &b_geninfo_mcweight);
   fChain->SetBranchAddress("NUP", &NUP, &b_NUP);
   fChain->SetBranchAddress("geninfo_htgen", &geninfo_htgen, &b_geninfo_htgen);
   fChain->SetBranchAddress("geninvmass", &geninvmass, &b_geninvmass);
   fChain->SetBranchAddress("weight_gen", &weight_gen, &b_weight_gen);
   fChain->SetBranchAddress("genmet_pt", &genmet_pt, &b_genmet_pt);
   fChain->SetBranchAddress("genmet_px", &genmet_px, &b_genmet_px);
   fChain->SetBranchAddress("genmet_py", &genmet_py, &b_genmet_py);
   fChain->SetBranchAddress("genmet_phi", &genmet_phi, &b_genmet_phi);
   fChain->SetBranchAddress("vbf_mjj", &vbf_mjj, &b_vbf_mjj);
   fChain->SetBranchAddress("vbf_deta", &vbf_deta, &b_vbf_deta);
   fChain->SetBranchAddress("vbf_n_central20", &vbf_n_central20, &b_vbf_n_central20);
   fChain->SetBranchAddress("vbf_n_central", &vbf_n_central, &b_vbf_n_central);
   fChain->SetBranchAddress("vbf_jdphi", &vbf_jdphi, &b_vbf_jdphi);
   fChain->SetBranchAddress("vbf_dijetpt", &vbf_dijetpt, &b_vbf_dijetpt);
   fChain->SetBranchAddress("vbf_dijetphi", &vbf_dijetphi, &b_vbf_dijetphi);
   fChain->SetBranchAddress("vbf_dphidijethiggs", &vbf_dphidijethiggs, &b_vbf_dphidijethiggs);
   fChain->SetBranchAddress("vbf_mindetajetvis", &vbf_mindetajetvis, &b_vbf_mindetajetvis);
   fChain->SetBranchAddress("jpt_1", &jpt_1, &b_jpt_1);
   fChain->SetBranchAddress("jeta_1", &jeta_1, &b_jeta_1);
   fChain->SetBranchAddress("jphi_1", &jphi_1, &b_jphi_1);
   fChain->SetBranchAddress("jet1_charge", &jet1_charge, &b_jet1_charge);
   fChain->SetBranchAddress("jet1_mass", &jet1_mass, &b_jet1_mass);
   fChain->SetBranchAddress("jmva_1", &jmva_1, &b_jmva_1);
   fChain->SetBranchAddress("jpuid_1", &jpuid_1, &b_jpuid_1);
   fChain->SetBranchAddress("jet1_flavour_parton", &jet1_flavour_parton, &b_jet1_flavour_parton);
   fChain->SetBranchAddress("jcsv_1", &jcsv_1, &b_jcsv_1);
   fChain->SetBranchAddress("jet1_genjet_pt", &jet1_genjet_pt, &b_jet1_genjet_pt);
   fChain->SetBranchAddress("jpt_2", &jpt_2, &b_jpt_2);
   fChain->SetBranchAddress("jeta_2", &jeta_2, &b_jeta_2);
   fChain->SetBranchAddress("jphi_2", &jphi_2, &b_jphi_2);
   fChain->SetBranchAddress("jet2_charge", &jet2_charge, &b_jet2_charge);
   fChain->SetBranchAddress("jet2_mass", &jet2_mass, &b_jet2_mass);
   fChain->SetBranchAddress("jmva_2", &jmva_2, &b_jmva_2);
   fChain->SetBranchAddress("jpuid_2", &jpuid_2, &b_jpuid_2);
   fChain->SetBranchAddress("jet2_flavour_parton", &jet2_flavour_parton, &b_jet2_flavour_parton);
   fChain->SetBranchAddress("jcsv_2", &jcsv_2, &b_jcsv_2);
   fChain->SetBranchAddress("jet2_genjet_pt", &jet2_genjet_pt, &b_jet2_genjet_pt);
   fChain->SetBranchAddress("bpt_1", &bpt_1, &b_bpt_1);
   fChain->SetBranchAddress("beta_1", &beta_1, &b_beta_1);
   fChain->SetBranchAddress("bphi_1", &bphi_1, &b_bphi_1);
   fChain->SetBranchAddress("bjet1_charge", &bjet1_charge, &b_bjet1_charge);
   fChain->SetBranchAddress("bjet1_mass", &bjet1_mass, &b_bjet1_mass);
   fChain->SetBranchAddress("bmva_1", &bmva_1, &b_bmva_1);
   fChain->SetBranchAddress("bpuid_1", &bpuid_1, &b_bpuid_1);
   fChain->SetBranchAddress("bjet1_flavour_parton", &bjet1_flavour_parton, &b_bjet1_flavour_parton);
   fChain->SetBranchAddress("bcsv_1", &bcsv_1, &b_bcsv_1);
   fChain->SetBranchAddress("bjet1_genjet_pt", &bjet1_genjet_pt, &b_bjet1_genjet_pt);
   fChain->SetBranchAddress("bpt_2", &bpt_2, &b_bpt_2);
   fChain->SetBranchAddress("beta_2", &beta_2, &b_beta_2);
   fChain->SetBranchAddress("bphi_2", &bphi_2, &b_bphi_2);
   fChain->SetBranchAddress("bjet2_charge", &bjet2_charge, &b_bjet2_charge);
   fChain->SetBranchAddress("bjet2_mass", &bjet2_mass, &b_bjet2_mass);
   fChain->SetBranchAddress("bmva_2", &bmva_2, &b_bmva_2);
   fChain->SetBranchAddress("bpuid_2", &bpuid_2, &b_bpuid_2);
   fChain->SetBranchAddress("bjet2_flavour_parton", &bjet2_flavour_parton, &b_bjet2_flavour_parton);
   fChain->SetBranchAddress("bcsv_2", &bcsv_2, &b_bcsv_2);
   fChain->SetBranchAddress("bjet2_genjet_pt", &bjet2_genjet_pt, &b_bjet2_genjet_pt);
   fChain->SetBranchAddress("HT_allJets", &HT_allJets, &b_HT_allJets);
   fChain->SetBranchAddress("HT_jets", &HT_jets, &b_HT_jets);
   fChain->SetBranchAddress("HT_bJets", &HT_bJets, &b_HT_bJets);
   fChain->SetBranchAddress("HT_cleanJets", &HT_cleanJets, &b_HT_cleanJets);
   fChain->SetBranchAddress("HT_jets30", &HT_jets30, &b_HT_jets30);
   fChain->SetBranchAddress("HT_cleanJets30", &HT_cleanJets30, &b_HT_cleanJets30);
   fChain->SetBranchAddress("genboson_pt", &genboson_pt, &b_genboson_pt);
   fChain->SetBranchAddress("genboson_eta", &genboson_eta, &b_genboson_eta);
   fChain->SetBranchAddress("genboson_phi", &genboson_phi, &b_genboson_phi);
   fChain->SetBranchAddress("genboson_charge", &genboson_charge, &b_genboson_charge);
   fChain->SetBranchAddress("genboson_mass", &genboson_mass, &b_genboson_mass);
   fChain->SetBranchAddress("genboson_pdgId", &genboson_pdgId, &b_genboson_pdgId);
   fChain->SetBranchAddress("gen_top_1_pt", &gen_top_1_pt, &b_gen_top_1_pt);
   fChain->SetBranchAddress("gen_top_2_pt", &gen_top_2_pt, &b_gen_top_2_pt);
   fChain->SetBranchAddress("gen_top_weight", &gen_top_weight, &b_gen_top_weight);
   fChain->SetBranchAddress("puppimet", &puppimet, &b_puppimet);
   fChain->SetBranchAddress("puppimetphi", &puppimetphi, &b_puppimetphi);
   fChain->SetBranchAddress("puppimt_1", &puppimt_1, &b_puppimt_1);
   fChain->SetBranchAddress("puppimt_2", &puppimt_2, &b_puppimt_2);
   fChain->SetBranchAddress("met", &met, &b_met);
   fChain->SetBranchAddress("metphi", &metphi, &b_metphi);
   fChain->SetBranchAddress("pfmet_mt1", &pfmet_mt1, &b_pfmet_mt1);
   fChain->SetBranchAddress("pfmet_mt2", &pfmet_mt2, &b_pfmet_mt2);
   fChain->SetBranchAddress("pt_2", &pt_2, &b_pt_2);
   fChain->SetBranchAddress("eta_2", &eta_2, &b_eta_2);
   fChain->SetBranchAddress("phi_2", &phi_2, &b_phi_2);
   fChain->SetBranchAddress("q_2", &q_2, &b_q_2);
   fChain->SetBranchAddress("m_2", &m_2, &b_m_2);
   fChain->SetBranchAddress("l2_jet_pt", &l2_jet_pt, &b_l2_jet_pt);
   fChain->SetBranchAddress("l2_jet_eta", &l2_jet_eta, &b_l2_jet_eta);
   fChain->SetBranchAddress("l2_jet_phi", &l2_jet_phi, &b_l2_jet_phi);
   fChain->SetBranchAddress("l2_jet_charge", &l2_jet_charge, &b_l2_jet_charge);
   fChain->SetBranchAddress("l2_jet_mass", &l2_jet_mass, &b_l2_jet_mass);
   fChain->SetBranchAddress("d0_2", &d0_2, &b_d0_2);
   fChain->SetBranchAddress("l2_dxy_error", &l2_dxy_error, &b_l2_dxy_error);
   fChain->SetBranchAddress("dZ_2", &dZ_2, &b_dZ_2);
   fChain->SetBranchAddress("l2_dz_error", &l2_dz_error, &b_l2_dz_error);
   fChain->SetBranchAddress("l2_weight", &l2_weight, &b_l2_weight);
   fChain->SetBranchAddress("trigweight_2", &trigweight_2, &b_trigweight_2);
   fChain->SetBranchAddress("l2_weight_eff_data_trigger", &l2_weight_eff_data_trigger, &b_l2_weight_eff_data_trigger);
   fChain->SetBranchAddress("l2_eff_trigger_data", &l2_eff_trigger_data, &b_l2_eff_trigger_data);
   fChain->SetBranchAddress("l2_eff_trigger_mc", &l2_eff_trigger_mc, &b_l2_eff_trigger_mc);
   fChain->SetBranchAddress("l2_weight_idiso", &l2_weight_idiso, &b_l2_weight_idiso);
   fChain->SetBranchAddress("l2_eff_idiso_data", &l2_eff_idiso_data, &b_l2_eff_idiso_data);
   fChain->SetBranchAddress("l2_eff_idiso_mc", &l2_eff_idiso_mc, &b_l2_eff_idiso_mc);
   fChain->SetBranchAddress("gen_match_2", &gen_match_2, &b_gen_match_2);
   fChain->SetBranchAddress("l2_decayMode", &l2_decayMode, &b_l2_decayMode);
   fChain->SetBranchAddress("l2_zImpact", &l2_zImpact, &b_l2_zImpact);
   fChain->SetBranchAddress("l2_dz_selfvertex", &l2_dz_selfvertex, &b_l2_dz_selfvertex);
   fChain->SetBranchAddress("l2_ptScale", &l2_ptScale, &b_l2_ptScale);
   fChain->SetBranchAddress("l2_againstElectronMVA6", &l2_againstElectronMVA6, &b_l2_againstElectronMVA6);
   fChain->SetBranchAddress("l2_againstMuon3", &l2_againstMuon3, &b_l2_againstMuon3);
   fChain->SetBranchAddress("byCombinedIsolationDeltaBetaCorrRaw3Hits_2", &byCombinedIsolationDeltaBetaCorrRaw3Hits_2, &b_byCombinedIsolationDeltaBetaCorrRaw3Hits_2);
   fChain->SetBranchAddress("l2_byIsolationMVArun2v1DBoldDMwLTraw", &l2_byIsolationMVArun2v1DBoldDMwLTraw, &b_l2_byIsolationMVArun2v1DBoldDMwLTraw);
   fChain->SetBranchAddress("l2_byIsolationMVArun2v1DBnewDMwLTraw", &l2_byIsolationMVArun2v1DBnewDMwLTraw, &b_l2_byIsolationMVArun2v1DBnewDMwLTraw);
   fChain->SetBranchAddress("l2_byIsolationMVArun2v1DBdR03oldDMwLTraw", &l2_byIsolationMVArun2v1DBdR03oldDMwLTraw, &b_l2_byIsolationMVArun2v1DBdR03oldDMwLTraw);
   fChain->SetBranchAddress("l2_byCombinedIsolationDeltaBetaCorr3Hits", &l2_byCombinedIsolationDeltaBetaCorr3Hits, &b_l2_byCombinedIsolationDeltaBetaCorr3Hits);
   fChain->SetBranchAddress("l2_byIsolationMVArun2v1DBoldDMwLT", &l2_byIsolationMVArun2v1DBoldDMwLT, &b_l2_byIsolationMVArun2v1DBoldDMwLT);
   fChain->SetBranchAddress("l2_byIsolationMVArun2v1DBnewDMwLT", &l2_byIsolationMVArun2v1DBnewDMwLT, &b_l2_byIsolationMVArun2v1DBnewDMwLT);
   fChain->SetBranchAddress("l2_byIsolationMVArun2v1DBdR03oldDMwLT", &l2_byIsolationMVArun2v1DBdR03oldDMwLT, &b_l2_byIsolationMVArun2v1DBdR03oldDMwLT);
   fChain->SetBranchAddress("chargedIsoPtSum_2", &chargedIsoPtSum_2, &b_chargedIsoPtSum_2);
   fChain->SetBranchAddress("decayModeFindingOldDMs_2", &decayModeFindingOldDMs_2, &b_decayModeFindingOldDMs_2);
   fChain->SetBranchAddress("l2_footprintCorrection", &l2_footprintCorrection, &b_l2_footprintCorrection);
   fChain->SetBranchAddress("neutralIsoPtSum_2", &neutralIsoPtSum_2, &b_neutralIsoPtSum_2);
   fChain->SetBranchAddress("puCorrPtSum_2", &puCorrPtSum_2, &b_puCorrPtSum_2);
   fChain->SetBranchAddress("l2_photonPtSumOutsideSignalCone", &l2_photonPtSumOutsideSignalCone, &b_l2_photonPtSumOutsideSignalCone);
   fChain->SetBranchAddress("mva_olddm_tight_2", &mva_olddm_tight_2, &b_mva_olddm_tight_2);
   fChain->SetBranchAddress("pt_1", &pt_1, &b_pt_1);
   fChain->SetBranchAddress("eta_1", &eta_1, &b_eta_1);
   fChain->SetBranchAddress("phi_1", &phi_1, &b_phi_1);
   fChain->SetBranchAddress("q_1", &q_1, &b_q_1);
   fChain->SetBranchAddress("m_1", &m_1, &b_m_1);
   fChain->SetBranchAddress("l1_jet_pt", &l1_jet_pt, &b_l1_jet_pt);
   fChain->SetBranchAddress("l1_jet_eta", &l1_jet_eta, &b_l1_jet_eta);
   fChain->SetBranchAddress("l1_jet_phi", &l1_jet_phi, &b_l1_jet_phi);
   fChain->SetBranchAddress("l1_jet_charge", &l1_jet_charge, &b_l1_jet_charge);
   fChain->SetBranchAddress("l1_jet_mass", &l1_jet_mass, &b_l1_jet_mass);
   fChain->SetBranchAddress("d0_1", &d0_1, &b_d0_1);
   fChain->SetBranchAddress("l1_dxy_error", &l1_dxy_error, &b_l1_dxy_error);
   fChain->SetBranchAddress("dZ_1", &dZ_1, &b_dZ_1);
   fChain->SetBranchAddress("l1_dz_error", &l1_dz_error, &b_l1_dz_error);
   fChain->SetBranchAddress("l1_weight", &l1_weight, &b_l1_weight);
   fChain->SetBranchAddress("trigweight_1", &trigweight_1, &b_trigweight_1);
   fChain->SetBranchAddress("l1_weight_eff_data_trigger", &l1_weight_eff_data_trigger, &b_l1_weight_eff_data_trigger);
   fChain->SetBranchAddress("l1_eff_trigger_data", &l1_eff_trigger_data, &b_l1_eff_trigger_data);
   fChain->SetBranchAddress("l1_eff_trigger_mc", &l1_eff_trigger_mc, &b_l1_eff_trigger_mc);
   fChain->SetBranchAddress("idisoweight_1", &idisoweight_1, &b_idisoweight_1);
   fChain->SetBranchAddress("l1_eff_idiso_data", &l1_eff_idiso_data, &b_l1_eff_idiso_data);
   fChain->SetBranchAddress("l1_eff_idiso_mc", &l1_eff_idiso_mc, &b_l1_eff_idiso_mc);
   fChain->SetBranchAddress("gen_match_1", &gen_match_1, &b_gen_match_1);
   fChain->SetBranchAddress("iso_1", &iso_1, &b_iso_1);
   //fChain->SetBranchAddress("l1_reliso05_04", &l1_reliso05_04, &b_l1_reliso05_04);
   fChain->SetBranchAddress("id_m_loose_1", &id_m_loose_1, &b_id_m_loose_1);
   fChain->SetBranchAddress("id_m_medium_1", &id_m_medium_1, &b_id_m_medium_1);
   fChain->SetBranchAddress("id_m_tight_1", &id_m_tight_1, &b_id_m_tight_1);
   fChain->SetBranchAddress("id_m_tightnovtx_1", &id_m_tightnovtx_1, &b_id_m_tightnovtx_1);
   fChain->SetBranchAddress("id_m_highpt_1", &id_m_highpt_1, &b_id_m_highpt_1);
   fChain->SetBranchAddress("l1_dxy_innertrack", &l1_dxy_innertrack, &b_l1_dxy_innertrack);
   fChain->SetBranchAddress("l1_dz_innertrack", &l1_dz_innertrack, &b_l1_dz_innertrack);
   fChain->SetBranchAddress("l2_gen_pt", &l2_gen_pt, &b_l2_gen_pt);
   fChain->SetBranchAddress("l2_gen_eta", &l2_gen_eta, &b_l2_gen_eta);
   fChain->SetBranchAddress("l2_gen_phi", &l2_gen_phi, &b_l2_gen_phi);
   fChain->SetBranchAddress("l2_gen_charge", &l2_gen_charge, &b_l2_gen_charge);
   fChain->SetBranchAddress("l2_gen_mass", &l2_gen_mass, &b_l2_gen_mass);
   fChain->SetBranchAddress("l2_gen_pdgId", &l2_gen_pdgId, &b_l2_gen_pdgId);
   fChain->SetBranchAddress("l2_gen_lepfromtau", &l2_gen_lepfromtau, &b_l2_gen_lepfromtau);
   fChain->SetBranchAddress("l1_gen_pt", &l1_gen_pt, &b_l1_gen_pt);
   fChain->SetBranchAddress("l1_gen_eta", &l1_gen_eta, &b_l1_gen_eta);
   fChain->SetBranchAddress("l1_gen_phi", &l1_gen_phi, &b_l1_gen_phi);
   fChain->SetBranchAddress("l1_gen_charge", &l1_gen_charge, &b_l1_gen_charge);
   fChain->SetBranchAddress("l1_gen_mass", &l1_gen_mass, &b_l1_gen_mass);
   fChain->SetBranchAddress("l1_gen_pdgId", &l1_gen_pdgId, &b_l1_gen_pdgId);
   fChain->SetBranchAddress("l1_gen_lepfromtau", &l1_gen_lepfromtau, &b_l1_gen_lepfromtau);
   fChain->SetBranchAddress("l2_gen_vis_pt", &l2_gen_vis_pt, &b_l2_gen_vis_pt);
   fChain->SetBranchAddress("l2_gen_vis_eta", &l2_gen_vis_eta, &b_l2_gen_vis_eta);
   fChain->SetBranchAddress("l2_gen_vis_phi", &l2_gen_vis_phi, &b_l2_gen_vis_phi);
   fChain->SetBranchAddress("l2_gen_vis_charge", &l2_gen_vis_charge, &b_l2_gen_vis_charge);
   fChain->SetBranchAddress("l2_gen_vis_mass", &l2_gen_vis_mass, &b_l2_gen_vis_mass);
   fChain->SetBranchAddress("l2_gen_decaymode", &l2_gen_decaymode, &b_l2_gen_decaymode);
   fChain->SetBranchAddress("l2_gen_nc_ratio", &l2_gen_nc_ratio, &b_l2_gen_nc_ratio);
   fChain->SetBranchAddress("l2_nc_ratio", &l2_nc_ratio, &b_l2_nc_ratio);
   fChain->SetBranchAddress("l2_weight_fakerate", &l2_weight_fakerate, &b_l2_weight_fakerate);
   fChain->SetBranchAddress("l2_weight_fakerate_up", &l2_weight_fakerate_up, &b_l2_weight_fakerate_up);
   fChain->SetBranchAddress("l2_weight_fakerate_down", &l2_weight_fakerate_down, &b_l2_weight_fakerate_down);
   
   //pileup
   hNPV = new TH1F(("hNPV_"+suffix_).c_str(),"",45,0.0,45.0);
   hNPU = new TH1F(("hNPU_"+suffix_).c_str(),"",40,0.0,40.0);
   //muon
   hMuPt = new TH1F(("hMuPt_"+suffix_).c_str(),"",100,0.0,750.0);
   hMuPhi = new TH1F(("hMuPhi_"+suffix_).c_str(),"",100,-TMath::Pi(),TMath::Pi());
   hMuEta = new TH1F(("hMuEta_"+suffix_).c_str(),"",100,-2.5,2.5);
   hMud0 = new TH1F(("hMud0_"+suffix_).c_str(),"",100,-0.05,0.05);
   hMudZ = new TH1F(("hMudZ_"+suffix_).c_str(),"",100,-0.2,0.2);
   hMuIso = new TH1F(("hMuIso_"+suffix_).c_str(),"",100,0,25);
   hMuIsoZoom = new TH1F(("hMuIsoZoom_"+suffix_).c_str(),"",100,0,0.5);
   hMuMt = new TH1F(("hMuMt_"+suffix_).c_str(),"",100,0,450);
   hMuPFMt = new TH1F(("hMuPFMt_"+suffix_).c_str(),"",100,0,450);
   hMuM = new TH1F(("hMuM_"+suffix_).c_str(),"",100,-0.6,0.4);
   hMuGenMatch = new TH1F(("hMuGenMatch_"+suffix_).c_str(),"",6,1,7);
   //tau
   hTauPt = new TH1F(("hTauPt_"+suffix_).c_str(),"",100,0,1300);
   hTauPhi = new TH1F(("hTauPhi_"+suffix_).c_str(),"",100,-TMath::Pi(),TMath::Pi());
   hTaud0 = new TH1F(("hTaud0_"+suffix_).c_str(),"",100,-0.25,0.25);
   hTaudZ = new TH1F(("hTaudZ_"+suffix_).c_str(),"",100,-0.2,0.2);
   hTauIso = new TH1F(("hTauIso_"+suffix_).c_str(),"",100,0,10);
   hTauMt = new TH1F(("hTauMt_"+suffix_).c_str(),"",100,0,1800);
   hTauPFMt = new TH1F(("hTauPFMt_"+suffix_).c_str(),"",100,0,1800);
   hTauEta = new TH1F(("hTauEta_"+suffix_).c_str(),"",100,-2.5,2.5);
   hTauGenMatch = new TH1F(("hTauGenMatch_"+suffix_).c_str(),"",6,1,7);
   //jets
   hJet1Pt = new TH1F(("hLeadingJetPt_"+suffix_).c_str(),"",100,0,1500);
   hJet1Phi = new TH1F(("hLeadingJetPhi_"+suffix_).c_str(),"",100,-TMath::Pi(),TMath::Pi());
   hJet1Eta = new TH1F(("hLeadingJetEta_"+suffix_).c_str(),"",100,-5,5);
   hJet2Pt = new TH1F(("hTrailingJetPt_"+suffix_).c_str(),"",100,0,1100);
   hJet2Phi = new TH1F(("hTrailingJetPhi_"+suffix_).c_str(),"",100,-TMath::Pi(),TMath::Pi());
   hJet2Eta = new TH1F(("hTrailingJetEta_"+suffix_).c_str(),"",100,-5,5);
   hBJet1Pt = new TH1F(("hLeadingBJetPt_"+suffix_).c_str(),"",100,0,1500);
   hBJet1Phi = new TH1F(("hLeadingBJetPhi_"+suffix_).c_str(),"",100,-TMath::Pi(),TMath::Pi());
   hBJet1Eta = new TH1F(("hLeadingBJetEta_"+suffix_).c_str(),"",100,-5,5);
   hBJet2Pt = new TH1F(("hTrailingBJetPt_"+suffix_).c_str(),"",100,0,1100);
   hBJet2Phi = new TH1F(("hTrailingBJetPhi_"+suffix_).c_str(),"",100,-TMath::Pi(),TMath::Pi());
   hBJet2Eta = new TH1F(("hTrailingBJetEta_"+suffix_).c_str(),"",100,-5,5);
   //MET
   hMET = new TH1F(("hMET_"+suffix_).c_str(),"",100,0,1300);
   hMETphi = new TH1F(("hMETphi_"+suffix_).c_str(),"",100,-TMath::Pi(),TMath::Pi());
   hMVAMET = new TH1F(("hMVAMET_"+suffix_).c_str(),"",100,0,1300);
   hMVAMETphi = new TH1F(("hMVAMETphi_"+suffix_).c_str(),"",100,-TMath::Pi(),TMath::Pi());
   hMVACov00 = new TH1F(("hMVACov00_"+suffix_).c_str(),"",100,0,17000);
   hMVACov01 = new TH1F(("hMVACov01_"+suffix_).c_str(),"",100,-5000,6000);
   hMVACov10 = new TH1F(("hMVACov10_"+suffix_).c_str(),"",100,-5000,6000);
   hMVACov11 = new TH1F(("hMVACov11_"+suffix_).c_str(),"",100,0,18000);
   hPZetaVis = new TH1F(("hPZetaVis_"+suffix_).c_str(),"",100,0,750);
   hPZetaMiss = new TH1F(("hPZetaMiss_"+suffix_).c_str(),"",100,-350,750);
   hPFPZetaMiss = new TH1F(("hPFPZetaMiss_"+suffix_).c_str(),"",100,-400,800);
   //di-tau
   hPtTT = new TH1F(("hPtTT_"+suffix_).c_str(),"",100,0,1300);
   hMtTot = new TH1F(("hMtTot_"+suffix_).c_str(),"",100,0,1700);
   hMVis = new TH1F(("hMVis_"+suffix_).c_str(),"",100,0,1700);
   hMSV = new TH1F(("hMSV_"+suffix_).c_str(),"",100,-1,1);
   //vbf system
   hMJJ = new TH1F(("hMJJ_"+suffix_).c_str(),"",100,0,6000);
   hJdEta = new TH1F(("hJdEta_"+suffix_).c_str(),"",100,0,10);
   hJdPhi = new TH1F(("hJdPhi_"+suffix_).c_str(),"",100,0,TMath::Pi());
   //extra lepton vetos
   hDiLeptonVeto = new TH1F(("hDiLeptonVeto_"+suffix_).c_str(),"",2,0,2);
   hExtraMuonVeto = new TH1F(("hExtraMuonVeto_"+suffix_).c_str(),"",2,0,2);
   hExtraElectronVeto = new TH1F(("hExtraElectronVeto_"+suffix_).c_str(),"",2,0,2);
   
   Notify();
}

Bool_t SynchCERN::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void SynchCERN::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t SynchCERN::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef SynchCERN_cxx
