//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Oct  5 12:51:45 2016 by ROOT version 6.06/00
// from TTree mt/mutau
// found on file: /scratch/apyskir/Data/Synch/VBFHToTauTauM125.root
//////////////////////////////////////////////////////////

#ifndef SynchLIP_h
#define SynchLIP_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class SynchLIP {
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
   TH1F*hTauGenMatch;
   //Jets
   TH1F*hJet1Pt;
   TH1F*hJet1Phi;
   TH1F*hJet1Eta;
   TH1F*hJet2Pt;
   TH1F*hJet2Phi;
   TH1F*hJet2Eta;
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

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           NUP;
   Int_t           dielectron_veto;
   Int_t           dimuon_veto;
   Int_t           evt;
   Int_t           extraelec_veto;
   Int_t           extramuon_veto;
   Int_t           gen_match_1;
   Int_t           gen_match_2;
   Int_t           lumi;
   Int_t           njets;
   Int_t           njetspt20;
   Int_t           npu;
   Int_t           npv;
   Int_t           run;
   Double_t        abs_mt_tot;
   Double_t        againstElectronLooseMVA6_1;
   Double_t        againstElectronLooseMVA6_2;
   Double_t        againstElectronMediumMVA6_1;
   Double_t        againstElectronMediumMVA6_2;
   Double_t        againstElectronTightMVA6_1;
   Double_t        againstElectronTightMVA6_2;
   Double_t        againstElectronVLooseMVA6_1;
   Double_t        againstElectronVLooseMVA6_2;
   Double_t        againstElectronVTightMVA6_1;
   Double_t        againstElectronVTightMVA6_2;
   Double_t        againstMuonLoose3_1;
   Double_t        againstMuonLoose3_2;
   Double_t        againstMuonTight3_1;
   Double_t        againstMuonTight3_2;
   Double_t        byCombinedIsolationDeltaBetaCorrRaw3Hits_1;
   Double_t        byCombinedIsolationDeltaBetaCorrRaw3Hits_2;
   Double_t        byIsolationMVA3newDMwLTraw_1;
   Double_t        byIsolationMVA3newDMwLTraw_2;
   Double_t        byIsolationMVA3newDMwoLTraw_1;
   Double_t        byIsolationMVA3newDMwoLTraw_2;
   Double_t        byIsolationMVA3oldDMwLTraw_1;
   Double_t        byIsolationMVA3oldDMwLTraw_2;
   Double_t        byIsolationMVA3oldDMwoLTraw_1;
   Double_t        byIsolationMVA3oldDMwoLTraw_2;
   Double_t        byLooseCombinedIsolationDeltaBetaCorr3Hits_1;
   Double_t        byLooseCombinedIsolationDeltaBetaCorr3Hits_2;
   Double_t        byLooseIsolationMVArun2v1DBoldDMwLT_1;
   Double_t        byLooseIsolationMVArun2v1DBoldDMwLT_2;
   Double_t        byMediumCombinedIsolationDeltaBetaCorr3Hits_1;
   Double_t        byMediumCombinedIsolationDeltaBetaCorr3Hits_2;
   Double_t        byMediumIsolationMVArun2v1DBoldDMwLT_1;
   Double_t        byMediumIsolationMVArun2v1DBoldDMwLT_2;
   Double_t        byTightCombinedIsolationDeltaBetaCorr3Hits_1;
   Double_t        byTightCombinedIsolationDeltaBetaCorr3Hits_2;
   Double_t        byTightIsolationMVArun2v1DBoldDMwLT_1;
   Double_t        byTightIsolationMVArun2v1DBoldDMwLT_2;
   Double_t        byVLooseIsolationMVArun2v1DBoldDMwLT_1;
   Double_t        byVLooseIsolationMVArun2v1DBoldDMwLT_2;
   Double_t        byVTightIsolationMVArun2v1DBoldDMwLT_1;
   Double_t        byVTightIsolationMVArun2v1DBoldDMwLT_2;
   Double_t        chargedIsoPtSum_1;
   Double_t        chargedIsoPtSum_2;
   Double_t        d0_1;
   Double_t        d0_2;
   Double_t        dZ_1;
   Double_t        dZ_2;
   Double_t        decayModeFindingOldDMs_1;
   Double_t        decayModeFindingOldDMs_2;
   Double_t        eta_1;
   Double_t        eta_2;
   Double_t        eta_sv;
   Double_t        id_e_cut_loose_1;
   Double_t        id_e_cut_loose_2;
   Double_t        id_e_cut_medium_1;
   Double_t        id_e_cut_medium_2;
   Double_t        id_e_cut_tight_1;
   Double_t        id_e_cut_tight_2;
   Double_t        id_e_cut_veto_1;
   Double_t        id_e_cut_veto_2;
   Double_t        id_e_mva_nt_loose_1;
   Double_t        id_m_highpt_1;
   Double_t        id_m_highpt_2;
   Double_t        id_m_loose_1;
   Double_t        id_m_loose_2;
   Double_t        id_m_medium_1;
   Double_t        id_m_medium_2;
   Double_t        id_m_tight_1;
   Double_t        id_m_tight_2;
   Double_t        id_m_tightnovtx_1;
   Double_t        id_m_tightnovtx_2;
   Double_t        idisoweight_1;
   Double_t        idisoweight_2;
   Double_t        iso_1;
   Double_t        iso_2;
   Double_t        jeta_1;
   Double_t        jmva_1;
   Double_t        jphi_1;
   Double_t        jpt_1;
   Double_t        jrawf_1;
   Double_t        m_1;
   Double_t        m_2;
   Double_t        m_sv;
   Double_t        m_vis;
   Double_t        met;
   Double_t        met_sv;
   Double_t        metcov00;
   Double_t        metcov01;
   Double_t        metcov10;
   Double_t        metcov11;
   Double_t        metphi;
   Double_t        mt_1;
   Double_t        mt_2;
   Double_t        mt_sv;
   Double_t        mvacov00;
   Double_t        mvacov01;
   Double_t        mvacov10;
   Double_t        mvacov11;
   Double_t        mvamet;
   Double_t        mvametphi;
   Double_t        neutralIsoPtSum_1;
   Double_t        neutralIsoPtSum_2;
   Double_t        pfmt_1;
   Double_t        pfpzetamiss;
   Double_t        phi_1;
   Double_t        phi_2;
   Double_t        phi_sv;
   Double_t        pt_1;
   Double_t        pt_2;
   Double_t        pt_sv;
   Double_t        pt_tt;
   Double_t        puCorrPtSum_1;
   Double_t        puCorrPtSum_2;
   Double_t        puppimet;
   Double_t        puppimetphi;
   Double_t        puppimt_1;
   Double_t        puweight;
   Double_t        pzetamiss;
   Double_t        pzetavis;
   Double_t        q_1;
   Double_t        q_2;
   Double_t        rho;
   Double_t        trigweight_1;
   Double_t        weight;

   // List of branches
   TBranch        *b_NUP;   //!
   TBranch        *b_dielectron_veto;   //!
   TBranch        *b_dimuon_veto;   //!
   TBranch        *b_evt;   //!
   TBranch        *b_extraelec_veto;   //!
   TBranch        *b_extramuon_veto;   //!
   TBranch        *b_gen_match_1;   //!
   TBranch        *b_gen_match_2;   //!
   TBranch        *b_lumi;   //!
   TBranch        *b_njets;   //!
   TBranch        *b_njetspt20;   //!
   TBranch        *b_npu;   //!
   TBranch        *b_npv;   //!
   TBranch        *b_run;   //!
   TBranch        *b_abs_mt_tot;   //!
   TBranch        *b_againstElectronLooseMVA6_1;   //!
   TBranch        *b_againstElectronLooseMVA6_2;   //!
   TBranch        *b_againstElectronMediumMVA6_1;   //!
   TBranch        *b_againstElectronMediumMVA6_2;   //!
   TBranch        *b_againstElectronTightMVA6_1;   //!
   TBranch        *b_againstElectronTightMVA6_2;   //!
   TBranch        *b_againstElectronVLooseMVA6_1;   //!
   TBranch        *b_againstElectronVLooseMVA6_2;   //!
   TBranch        *b_againstElectronVTightMVA6_1;   //!
   TBranch        *b_againstElectronVTightMVA6_2;   //!
   TBranch        *b_againstMuonLoose3_1;   //!
   TBranch        *b_againstMuonLoose3_2;   //!
   TBranch        *b_againstMuonTight3_1;   //!
   TBranch        *b_againstMuonTight3_2;   //!
   TBranch        *b_byCombinedIsolationDeltaBetaCorrRaw3Hits_1;   //!
   TBranch        *b_byCombinedIsolationDeltaBetaCorrRaw3Hits_2;   //!
   TBranch        *b_byIsolationMVA3newDMwLTraw_1;   //!
   TBranch        *b_byIsolationMVA3newDMwLTraw_2;   //!
   TBranch        *b_byIsolationMVA3newDMwoLTraw_1;   //!
   TBranch        *b_byIsolationMVA3newDMwoLTraw_2;   //!
   TBranch        *b_byIsolationMVA3oldDMwLTraw_1;   //!
   TBranch        *b_byIsolationMVA3oldDMwLTraw_2;   //!
   TBranch        *b_byIsolationMVA3oldDMwoLTraw_1;   //!
   TBranch        *b_byIsolationMVA3oldDMwoLTraw_2;   //!
   TBranch        *b_byLooseCombinedIsolationDeltaBetaCorr3Hits_1;   //!
   TBranch        *b_byLooseCombinedIsolationDeltaBetaCorr3Hits_2;   //!
   TBranch        *b_byLooseIsolationMVArun2v1DBoldDMwLT_1;   //!
   TBranch        *b_byLooseIsolationMVArun2v1DBoldDMwLT_2;   //!
   TBranch        *b_byMediumCombinedIsolationDeltaBetaCorr3Hits_1;   //!
   TBranch        *b_byMediumCombinedIsolationDeltaBetaCorr3Hits_2;   //!
   TBranch        *b_byMediumIsolationMVArun2v1DBoldDMwLT_1;   //!
   TBranch        *b_byMediumIsolationMVArun2v1DBoldDMwLT_2;   //!
   TBranch        *b_byTightCombinedIsolationDeltaBetaCorr3Hits_1;   //!
   TBranch        *b_byTightCombinedIsolationDeltaBetaCorr3Hits_2;   //!
   TBranch        *b_byTightIsolationMVArun2v1DBoldDMwLT_1;   //!
   TBranch        *b_byTightIsolationMVArun2v1DBoldDMwLT_2;   //!
   TBranch        *b_byVLooseIsolationMVArun2v1DBoldDMwLT_1;   //!
   TBranch        *b_byVLooseIsolationMVArun2v1DBoldDMwLT_2;   //!
   TBranch        *b_byVTightIsolationMVArun2v1DBoldDMwLT_1;   //!
   TBranch        *b_byVTightIsolationMVArun2v1DBoldDMwLT_2;   //!
   TBranch        *b_chargedIsoPtSum_1;   //!
   TBranch        *b_chargedIsoPtSum_2;   //!
   TBranch        *b_d0_1;   //!
   TBranch        *b_d0_2;   //!
   TBranch        *b_dZ_1;   //!
   TBranch        *b_dZ_2;   //!
   TBranch        *b_decayModeFindingOldDMs_1;   //!
   TBranch        *b_decayModeFindingOldDMs_2;   //!
   TBranch        *b_eta_1;   //!
   TBranch        *b_eta_2;   //!
   TBranch        *b_eta_sv;   //!
   TBranch        *b_id_e_cut_loose_1;   //!
   TBranch        *b_id_e_cut_loose_2;   //!
   TBranch        *b_id_e_cut_medium_1;   //!
   TBranch        *b_id_e_cut_medium_2;   //!
   TBranch        *b_id_e_cut_tight_1;   //!
   TBranch        *b_id_e_cut_tight_2;   //!
   TBranch        *b_id_e_cut_veto_1;   //!
   TBranch        *b_id_e_cut_veto_2;   //!
   TBranch        *b_id_e_mva_nt_loose_1;   //!
   TBranch        *b_id_m_highpt_1;   //!
   TBranch        *b_id_m_highpt_2;   //!
   TBranch        *b_id_m_loose_1;   //!
   TBranch        *b_id_m_loose_2;   //!
   TBranch        *b_id_m_medium_1;   //!
   TBranch        *b_id_m_medium_2;   //!
   TBranch        *b_id_m_tight_1;   //!
   TBranch        *b_id_m_tight_2;   //!
   TBranch        *b_id_m_tightnovtx_1;   //!
   TBranch        *b_id_m_tightnovtx_2;   //!
   TBranch        *b_idisoweight_1;   //!
   TBranch        *b_idisoweight_2;   //!
   TBranch        *b_iso_1;   //!
   TBranch        *b_iso_2;   //!
   TBranch        *b_jeta_1;   //!
   TBranch        *b_jmva_1;   //!
   TBranch        *b_jphi_1;   //!
   TBranch        *b_jpt_1;   //!
   TBranch        *b_jrawf_1;   //!
   TBranch        *b_m_1;   //!
   TBranch        *b_m_2;   //!
   TBranch        *b_m_sv;   //!
   TBranch        *b_m_vis;   //!
   TBranch        *b_met;   //!
   TBranch        *b_met_sv;   //!
   TBranch        *b_metcov00;   //!
   TBranch        *b_metcov01;   //!
   TBranch        *b_metcov10;   //!
   TBranch        *b_metcov11;   //!
   TBranch        *b_metphi;   //!
   TBranch        *b_mt_1;   //!
   TBranch        *b_mt_2;   //!
   TBranch        *b_mt_sv;   //!
   TBranch        *b_mvacov00;   //!
   TBranch        *b_mvacov01;   //!
   TBranch        *b_mvacov10;   //!
   TBranch        *b_mvacov11;   //!
   TBranch        *b_mvamet;   //!
   TBranch        *b_mvametphi;   //!
   TBranch        *b_neutralIsoPtSum_1;   //!
   TBranch        *b_neutralIsoPtSum_2;   //!
   TBranch        *b_pfmt_1;   //!
   TBranch        *b_pfpzetamiss;   //!
   TBranch        *b_phi_1;   //!
   TBranch        *b_phi_2;   //!
   TBranch        *b_phi_sv;   //!
   TBranch        *b_pt_1;   //!
   TBranch        *b_pt_2;   //!
   TBranch        *b_pt_sv;   //!
   TBranch        *b_pt_tt;   //!
   TBranch        *b_puCorrPtSum_1;   //!
   TBranch        *b_puCorrPtSum_2;   //!
   TBranch        *b_puppimet;   //!
   TBranch        *b_puppimetphi;   //!
   TBranch        *b_puppimt_1;   //!
   TBranch        *b_puweight;   //!
   TBranch        *b_pzetamiss;   //!
   TBranch        *b_pzetavis;   //!
   TBranch        *b_q_1;   //!
   TBranch        *b_q_2;   //!
   TBranch        *b_rho;   //!
   TBranch        *b_trigweight_1;   //!
   TBranch        *b_weight;   //!

   SynchLIP(TTree *tree=0, std::string suffix="");
   virtual ~SynchLIP();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef SynchLIP_cxx
SynchLIP::SynchLIP(TTree *tree, std::string suffix) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   suffix_ = suffix;
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/scratch/apyskir/Data/Synch/VBFHToTauTauM125.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("/scratch/apyskir/Data/Synch/VBFHToTauTauM125.root");
      }
      TDirectory * dir = (TDirectory*)f->Get("/scratch/apyskir/Data/Synch/VBFHToTauTauM125.root:/h2taus");
      dir->GetObject("mt",tree);

   }
   Init(tree);
}

SynchLIP::~SynchLIP()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t SynchLIP::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t SynchLIP::LoadTree(Long64_t entry)
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

void SynchLIP::Init(TTree *tree)
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

   fChain->SetBranchAddress("NUP", &NUP, &b_NUP);
   fChain->SetBranchAddress("dielectron_veto", &dielectron_veto, &b_dielectron_veto);
   fChain->SetBranchAddress("dimuon_veto", &dimuon_veto, &b_dimuon_veto);
   fChain->SetBranchAddress("evt", &evt, &b_evt);
   fChain->SetBranchAddress("extraelec_veto", &extraelec_veto, &b_extraelec_veto);
   fChain->SetBranchAddress("extramuon_veto", &extramuon_veto, &b_extramuon_veto);
   fChain->SetBranchAddress("gen_match_1", &gen_match_1, &b_gen_match_1);
   fChain->SetBranchAddress("gen_match_2", &gen_match_2, &b_gen_match_2);
   fChain->SetBranchAddress("lumi", &lumi, &b_lumi);
   fChain->SetBranchAddress("njets", &njets, &b_njets);
   fChain->SetBranchAddress("njetspt20", &njetspt20, &b_njetspt20);
   fChain->SetBranchAddress("npu", &npu, &b_npu);
   fChain->SetBranchAddress("npv", &npv, &b_npv);
   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("abs_mt_tot", &abs_mt_tot, &b_abs_mt_tot);
   fChain->SetBranchAddress("againstElectronLooseMVA6_1", &againstElectronLooseMVA6_1, &b_againstElectronLooseMVA6_1);
   fChain->SetBranchAddress("againstElectronLooseMVA6_2", &againstElectronLooseMVA6_2, &b_againstElectronLooseMVA6_2);
   fChain->SetBranchAddress("againstElectronMediumMVA6_1", &againstElectronMediumMVA6_1, &b_againstElectronMediumMVA6_1);
   fChain->SetBranchAddress("againstElectronMediumMVA6_2", &againstElectronMediumMVA6_2, &b_againstElectronMediumMVA6_2);
   fChain->SetBranchAddress("againstElectronTightMVA6_1", &againstElectronTightMVA6_1, &b_againstElectronTightMVA6_1);
   fChain->SetBranchAddress("againstElectronTightMVA6_2", &againstElectronTightMVA6_2, &b_againstElectronTightMVA6_2);
   fChain->SetBranchAddress("againstElectronVLooseMVA6_1", &againstElectronVLooseMVA6_1, &b_againstElectronVLooseMVA6_1);
   fChain->SetBranchAddress("againstElectronVLooseMVA6_2", &againstElectronVLooseMVA6_2, &b_againstElectronVLooseMVA6_2);
   fChain->SetBranchAddress("againstElectronVTightMVA6_1", &againstElectronVTightMVA6_1, &b_againstElectronVTightMVA6_1);
   fChain->SetBranchAddress("againstElectronVTightMVA6_2", &againstElectronVTightMVA6_2, &b_againstElectronVTightMVA6_2);
   fChain->SetBranchAddress("againstMuonLoose3_1", &againstMuonLoose3_1, &b_againstMuonLoose3_1);
   fChain->SetBranchAddress("againstMuonLoose3_2", &againstMuonLoose3_2, &b_againstMuonLoose3_2);
   fChain->SetBranchAddress("againstMuonTight3_1", &againstMuonTight3_1, &b_againstMuonTight3_1);
   fChain->SetBranchAddress("againstMuonTight3_2", &againstMuonTight3_2, &b_againstMuonTight3_2);
   fChain->SetBranchAddress("byCombinedIsolationDeltaBetaCorrRaw3Hits_1", &byCombinedIsolationDeltaBetaCorrRaw3Hits_1, &b_byCombinedIsolationDeltaBetaCorrRaw3Hits_1);
   fChain->SetBranchAddress("byCombinedIsolationDeltaBetaCorrRaw3Hits_2", &byCombinedIsolationDeltaBetaCorrRaw3Hits_2, &b_byCombinedIsolationDeltaBetaCorrRaw3Hits_2);
   fChain->SetBranchAddress("byIsolationMVA3newDMwLTraw_1", &byIsolationMVA3newDMwLTraw_1, &b_byIsolationMVA3newDMwLTraw_1);
   fChain->SetBranchAddress("byIsolationMVA3newDMwLTraw_2", &byIsolationMVA3newDMwLTraw_2, &b_byIsolationMVA3newDMwLTraw_2);
   fChain->SetBranchAddress("byIsolationMVA3newDMwoLTraw_1", &byIsolationMVA3newDMwoLTraw_1, &b_byIsolationMVA3newDMwoLTraw_1);
   fChain->SetBranchAddress("byIsolationMVA3newDMwoLTraw_2", &byIsolationMVA3newDMwoLTraw_2, &b_byIsolationMVA3newDMwoLTraw_2);
   fChain->SetBranchAddress("byIsolationMVA3oldDMwLTraw_1", &byIsolationMVA3oldDMwLTraw_1, &b_byIsolationMVA3oldDMwLTraw_1);
   fChain->SetBranchAddress("byIsolationMVA3oldDMwLTraw_2", &byIsolationMVA3oldDMwLTraw_2, &b_byIsolationMVA3oldDMwLTraw_2);
   fChain->SetBranchAddress("byIsolationMVA3oldDMwoLTraw_1", &byIsolationMVA3oldDMwoLTraw_1, &b_byIsolationMVA3oldDMwoLTraw_1);
   fChain->SetBranchAddress("byIsolationMVA3oldDMwoLTraw_2", &byIsolationMVA3oldDMwoLTraw_2, &b_byIsolationMVA3oldDMwoLTraw_2);
   fChain->SetBranchAddress("byLooseCombinedIsolationDeltaBetaCorr3Hits_1", &byLooseCombinedIsolationDeltaBetaCorr3Hits_1, &b_byLooseCombinedIsolationDeltaBetaCorr3Hits_1);
   fChain->SetBranchAddress("byLooseCombinedIsolationDeltaBetaCorr3Hits_2", &byLooseCombinedIsolationDeltaBetaCorr3Hits_2, &b_byLooseCombinedIsolationDeltaBetaCorr3Hits_2);
   fChain->SetBranchAddress("byLooseIsolationMVArun2v1DBoldDMwLT_1", &byLooseIsolationMVArun2v1DBoldDMwLT_1, &b_byLooseIsolationMVArun2v1DBoldDMwLT_1);
   fChain->SetBranchAddress("byLooseIsolationMVArun2v1DBoldDMwLT_2", &byLooseIsolationMVArun2v1DBoldDMwLT_2, &b_byLooseIsolationMVArun2v1DBoldDMwLT_2);
   fChain->SetBranchAddress("byMediumCombinedIsolationDeltaBetaCorr3Hits_1", &byMediumCombinedIsolationDeltaBetaCorr3Hits_1, &b_byMediumCombinedIsolationDeltaBetaCorr3Hits_1);
   fChain->SetBranchAddress("byMediumCombinedIsolationDeltaBetaCorr3Hits_2", &byMediumCombinedIsolationDeltaBetaCorr3Hits_2, &b_byMediumCombinedIsolationDeltaBetaCorr3Hits_2);
   fChain->SetBranchAddress("byMediumIsolationMVArun2v1DBoldDMwLT_1", &byMediumIsolationMVArun2v1DBoldDMwLT_1, &b_byMediumIsolationMVArun2v1DBoldDMwLT_1);
   fChain->SetBranchAddress("byMediumIsolationMVArun2v1DBoldDMwLT_2", &byMediumIsolationMVArun2v1DBoldDMwLT_2, &b_byMediumIsolationMVArun2v1DBoldDMwLT_2);
   fChain->SetBranchAddress("byTightCombinedIsolationDeltaBetaCorr3Hits_1", &byTightCombinedIsolationDeltaBetaCorr3Hits_1, &b_byTightCombinedIsolationDeltaBetaCorr3Hits_1);
   fChain->SetBranchAddress("byTightCombinedIsolationDeltaBetaCorr3Hits_2", &byTightCombinedIsolationDeltaBetaCorr3Hits_2, &b_byTightCombinedIsolationDeltaBetaCorr3Hits_2);
   fChain->SetBranchAddress("byTightIsolationMVArun2v1DBoldDMwLT_1", &byTightIsolationMVArun2v1DBoldDMwLT_1, &b_byTightIsolationMVArun2v1DBoldDMwLT_1);
   fChain->SetBranchAddress("byTightIsolationMVArun2v1DBoldDMwLT_2", &byTightIsolationMVArun2v1DBoldDMwLT_2, &b_byTightIsolationMVArun2v1DBoldDMwLT_2);
   fChain->SetBranchAddress("byVLooseIsolationMVArun2v1DBoldDMwLT_1", &byVLooseIsolationMVArun2v1DBoldDMwLT_1, &b_byVLooseIsolationMVArun2v1DBoldDMwLT_1);
   fChain->SetBranchAddress("byVLooseIsolationMVArun2v1DBoldDMwLT_2", &byVLooseIsolationMVArun2v1DBoldDMwLT_2, &b_byVLooseIsolationMVArun2v1DBoldDMwLT_2);
   fChain->SetBranchAddress("byVTightIsolationMVArun2v1DBoldDMwLT_1", &byVTightIsolationMVArun2v1DBoldDMwLT_1, &b_byVTightIsolationMVArun2v1DBoldDMwLT_1);
   fChain->SetBranchAddress("byVTightIsolationMVArun2v1DBoldDMwLT_2", &byVTightIsolationMVArun2v1DBoldDMwLT_2, &b_byVTightIsolationMVArun2v1DBoldDMwLT_2);
   fChain->SetBranchAddress("chargedIsoPtSum_1", &chargedIsoPtSum_1, &b_chargedIsoPtSum_1);
   fChain->SetBranchAddress("chargedIsoPtSum_2", &chargedIsoPtSum_2, &b_chargedIsoPtSum_2);
   fChain->SetBranchAddress("d0_1", &d0_1, &b_d0_1);
   fChain->SetBranchAddress("d0_2", &d0_2, &b_d0_2);
   fChain->SetBranchAddress("dZ_1", &dZ_1, &b_dZ_1);
   fChain->SetBranchAddress("dZ_2", &dZ_2, &b_dZ_2);
   fChain->SetBranchAddress("decayModeFindingOldDMs_1", &decayModeFindingOldDMs_1, &b_decayModeFindingOldDMs_1);
   fChain->SetBranchAddress("decayModeFindingOldDMs_2", &decayModeFindingOldDMs_2, &b_decayModeFindingOldDMs_2);
   fChain->SetBranchAddress("eta_1", &eta_1, &b_eta_1);
   fChain->SetBranchAddress("eta_2", &eta_2, &b_eta_2);
   fChain->SetBranchAddress("eta_sv", &eta_sv, &b_eta_sv);
   fChain->SetBranchAddress("id_e_cut_loose_1", &id_e_cut_loose_1, &b_id_e_cut_loose_1);
   fChain->SetBranchAddress("id_e_cut_loose_2", &id_e_cut_loose_2, &b_id_e_cut_loose_2);
   fChain->SetBranchAddress("id_e_cut_medium_1", &id_e_cut_medium_1, &b_id_e_cut_medium_1);
   fChain->SetBranchAddress("id_e_cut_medium_2", &id_e_cut_medium_2, &b_id_e_cut_medium_2);
   fChain->SetBranchAddress("id_e_cut_tight_1", &id_e_cut_tight_1, &b_id_e_cut_tight_1);
   fChain->SetBranchAddress("id_e_cut_tight_2", &id_e_cut_tight_2, &b_id_e_cut_tight_2);
   fChain->SetBranchAddress("id_e_cut_veto_1", &id_e_cut_veto_1, &b_id_e_cut_veto_1);
   fChain->SetBranchAddress("id_e_cut_veto_2", &id_e_cut_veto_2, &b_id_e_cut_veto_2);
   fChain->SetBranchAddress("id_e_mva_nt_loose_1", &id_e_mva_nt_loose_1, &b_id_e_mva_nt_loose_1);
   fChain->SetBranchAddress("id_m_highpt_1", &id_m_highpt_1, &b_id_m_highpt_1);
   fChain->SetBranchAddress("id_m_highpt_2", &id_m_highpt_2, &b_id_m_highpt_2);
   fChain->SetBranchAddress("id_m_loose_1", &id_m_loose_1, &b_id_m_loose_1);
   fChain->SetBranchAddress("id_m_loose_2", &id_m_loose_2, &b_id_m_loose_2);
   fChain->SetBranchAddress("id_m_medium_1", &id_m_medium_1, &b_id_m_medium_1);
   fChain->SetBranchAddress("id_m_medium_2", &id_m_medium_2, &b_id_m_medium_2);
   fChain->SetBranchAddress("id_m_tight_1", &id_m_tight_1, &b_id_m_tight_1);
   fChain->SetBranchAddress("id_m_tight_2", &id_m_tight_2, &b_id_m_tight_2);
   fChain->SetBranchAddress("id_m_tightnovtx_1", &id_m_tightnovtx_1, &b_id_m_tightnovtx_1);
   fChain->SetBranchAddress("id_m_tightnovtx_2", &id_m_tightnovtx_2, &b_id_m_tightnovtx_2);
   fChain->SetBranchAddress("idisoweight_1", &idisoweight_1, &b_idisoweight_1);
   fChain->SetBranchAddress("idisoweight_2", &idisoweight_2, &b_idisoweight_2);
   fChain->SetBranchAddress("iso_1", &iso_1, &b_iso_1);
   fChain->SetBranchAddress("iso_2", &iso_2, &b_iso_2);
   fChain->SetBranchAddress("jeta_1", &jeta_1, &b_jeta_1);
   fChain->SetBranchAddress("jmva_1", &jmva_1, &b_jmva_1);
   fChain->SetBranchAddress("jphi_1", &jphi_1, &b_jphi_1);
   fChain->SetBranchAddress("jpt_1", &jpt_1, &b_jpt_1);
   fChain->SetBranchAddress("jrawf_1", &jrawf_1, &b_jrawf_1);
   fChain->SetBranchAddress("m_1", &m_1, &b_m_1);
   fChain->SetBranchAddress("m_2", &m_2, &b_m_2);
   fChain->SetBranchAddress("m_sv", &m_sv, &b_m_sv);
   fChain->SetBranchAddress("m_vis", &m_vis, &b_m_vis);
   fChain->SetBranchAddress("met", &met, &b_met);
   fChain->SetBranchAddress("met_sv", &met_sv, &b_met_sv);
   fChain->SetBranchAddress("metcov00", &metcov00, &b_metcov00);
   fChain->SetBranchAddress("metcov01", &metcov01, &b_metcov01);
   fChain->SetBranchAddress("metcov10", &metcov10, &b_metcov10);
   fChain->SetBranchAddress("metcov11", &metcov11, &b_metcov11);
   fChain->SetBranchAddress("metphi", &metphi, &b_metphi);
   fChain->SetBranchAddress("mt_1", &mt_1, &b_mt_1);
   fChain->SetBranchAddress("mt_2", &mt_2, &b_mt_2);
   fChain->SetBranchAddress("mt_sv", &mt_sv, &b_mt_sv);
   fChain->SetBranchAddress("mvacov00", &mvacov00, &b_mvacov00);
   fChain->SetBranchAddress("mvacov01", &mvacov01, &b_mvacov01);
   fChain->SetBranchAddress("mvacov10", &mvacov10, &b_mvacov10);
   fChain->SetBranchAddress("mvacov11", &mvacov11, &b_mvacov11);
   fChain->SetBranchAddress("mvamet", &mvamet, &b_mvamet);
   fChain->SetBranchAddress("mvametphi", &mvametphi, &b_mvametphi);
   fChain->SetBranchAddress("neutralIsoPtSum_1", &neutralIsoPtSum_1, &b_neutralIsoPtSum_1);
   fChain->SetBranchAddress("neutralIsoPtSum_2", &neutralIsoPtSum_2, &b_neutralIsoPtSum_2);
   fChain->SetBranchAddress("pfmt_1", &pfmt_1, &b_pfmt_1);
   fChain->SetBranchAddress("pfpzetamiss", &pfpzetamiss, &b_pfpzetamiss);
   fChain->SetBranchAddress("phi_1", &phi_1, &b_phi_1);
   fChain->SetBranchAddress("phi_2", &phi_2, &b_phi_2);
   fChain->SetBranchAddress("phi_sv", &phi_sv, &b_phi_sv);
   fChain->SetBranchAddress("pt_1", &pt_1, &b_pt_1);
   fChain->SetBranchAddress("pt_2", &pt_2, &b_pt_2);
   fChain->SetBranchAddress("pt_sv", &pt_sv, &b_pt_sv);
   fChain->SetBranchAddress("pt_tt", &pt_tt, &b_pt_tt);
   fChain->SetBranchAddress("puCorrPtSum_1", &puCorrPtSum_1, &b_puCorrPtSum_1);
   fChain->SetBranchAddress("puCorrPtSum_2", &puCorrPtSum_2, &b_puCorrPtSum_2);
   fChain->SetBranchAddress("puppimet", &puppimet, &b_puppimet);
   fChain->SetBranchAddress("puppimetphi", &puppimetphi, &b_puppimetphi);
   fChain->SetBranchAddress("puppimt_1", &puppimt_1, &b_puppimt_1);
   fChain->SetBranchAddress("puweight", &puweight, &b_puweight);
   fChain->SetBranchAddress("pzetamiss", &pzetamiss, &b_pzetamiss);
   fChain->SetBranchAddress("pzetavis", &pzetavis, &b_pzetavis);
   fChain->SetBranchAddress("q_1", &q_1, &b_q_1);
   fChain->SetBranchAddress("q_2", &q_2, &b_q_2);
   fChain->SetBranchAddress("rho", &rho, &b_rho);
   fChain->SetBranchAddress("trigweight_1", &trigweight_1, &b_trigweight_1);
   fChain->SetBranchAddress("weight", &weight, &b_weight);
   
   //pileup
   hNPV = new TH1F(("hNPV_"+suffix_).c_str(),"",45,0.0,45.0);;
   hNPU = new TH1F(("hNPU_"+suffix_).c_str(),"",40,0.0,40.0);;;
   //muon
   hMuPt = new TH1F(("hMuPt_"+suffix_).c_str(),"",100,0.0,750.0);
   hMuPhi = new TH1F(("hMuPhi_"+suffix_).c_str(),"",100,-TMath::Pi(),TMath::Pi());
   hMuEta = new TH1F(("hMuEta_"+suffix_).c_str(),"",100,-2.5,2.5);
   hMud0 = new TH1F(("hMud0_"+suffix_).c_str(),"",100,-0.05,0.05);
   hMudZ = new TH1F(("hMudZ_"+suffix_).c_str(),"",100,-0.2,0.2);
   hMuIso = new TH1F(("hMuIso_"+suffix_).c_str(),"",100,0,25);
   hMuIsoZoom = new TH1F(("hMuIsoZoom_"+suffix_).c_str(),"",100,0,0.5);
   hMuMt = new TH1F(("hMuMt_"+suffix_).c_str(),"",100,0,450);
   hMuM = new TH1F(("hMuM_"+suffix_).c_str(),"",100,-0.6,0.4);
   hMuGenMatch = new TH1F(("hMuGenMatch_"+suffix_).c_str(),"",6,1,7);
   //tau
   hTauPt = new TH1F(("hTauPt_"+suffix_).c_str(),"",100,0,1300);
   hTauPhi = new TH1F(("hTauPhi_"+suffix_).c_str(),"",100,-TMath::Pi(),TMath::Pi());
   hTaud0 = new TH1F(("hTaud0_"+suffix_).c_str(),"",100,-0.25,0.25);
   hTaudZ = new TH1F(("hTaudZ_"+suffix_).c_str(),"",100,-0.2,0.2);
   hTauIso = new TH1F(("hTauIso_"+suffix_).c_str(),"",100,0,10);
   hTauMt = new TH1F(("hTauMt_"+suffix_).c_str(),"",100,0,1800);
   hTauEta = new TH1F(("hTauEta_"+suffix_).c_str(),"",100,-2.5,2.5);
   hTauGenMatch = new TH1F(("hTauGenMatch_"+suffix_).c_str(),"",6,1,7);
   //jets
   hJet1Pt = new TH1F(("hLeadingJetPt_"+suffix_).c_str(),"",100,0,1500);
   hJet1Phi = new TH1F(("hLeadingJetPhi_"+suffix_).c_str(),"",100,-TMath::Pi(),TMath::Pi());
   hJet1Eta = new TH1F(("hLeadingJetEta_"+suffix_).c_str(),"",100,-5,5);
   hJet2Pt = new TH1F(("hTrailingJetPt_"+suffix_).c_str(),"",100,0,1100);
   hJet2Phi = new TH1F(("hTrailingJetPhi_"+suffix_).c_str(),"",100,-TMath::Pi(),TMath::Pi());
   hJet2Eta = new TH1F(("hTrailingJetEta_"+suffix_).c_str(),"",100,-5,5);
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

Bool_t SynchLIP::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void SynchLIP::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t SynchLIP::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef SynchLIP_cxx
