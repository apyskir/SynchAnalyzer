//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Sep 13 13:25:33 2016 by ROOT version 6.06/01
// from TTree tree/Selections bit words
// found on file: RootAnalysis_SynchNTuple.root
//////////////////////////////////////////////////////////

#ifndef Synch_h
#define Synch_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include "TMath.h"

// Header file for the classes stored in the TTree if any.

class Synch {
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

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Float_t         eventWeight;
   ULong64_t       run;
   ULong64_t       lumi;
   ULong64_t       evt;
   ULong64_t       npv;
   ULong64_t       npu;
   Float_t         rho;
   Float_t         pt_1;
   Float_t         phi_1;
   Float_t         eta_1;
   Float_t         m_1;
   Float_t         q_1;
   Float_t         d0_1;
   Float_t         dZ_1;
   Float_t         mt_1;
   Float_t         pfmt_1;
   Float_t         puppimt_1;
   Float_t         iso_1;
   Float_t         id_e_mva_nt_loose_1;
   Int_t           gen_match_1;
   Float_t         againstElectronLooseMVA6_1;
   Float_t         againstElectronMediumMVA6_1;
   Float_t         againstElectronTightMVA6_1;
   Float_t         againstElectronVLooseMVA6_1;
   Float_t         againstElectronVTightMVA6_1;
   Float_t         againstMuonLoose3_1;
   Float_t         againstMuonTight3_1;
   Float_t         byCombinedIsolationDeltaBetaCorrRaw3Hits_1;
   Float_t         byIsolationMVA3newDMwoLTraw_1;
   Float_t         byIsolationMVA3oldDMwoLTraw_1;
   Float_t         byIsolationMVA3newDMwLTraw_1;
   Float_t         byIsolationMVA3oldDMwLTraw_1;
   Float_t         chargedIsoPtSum_1;
   Float_t         decayModeFindingOldDMs_1;
   Float_t         neutralIsoPtSum_1;
   Float_t         puCorrPtSum_1;
   Float_t         trigweight_1;
   Float_t         idisoweight_1;
   Float_t         pt_2;
   Float_t         phi_2;
   Float_t         eta_2;
   Float_t         m_2;
   Float_t         q_2;
   Float_t         d0_2;
   Float_t         dZ_2;
   Float_t         mt_2;
   Float_t         pfmt_2;
   Float_t         puppimt_2;
   Float_t         iso_2;
   Float_t         id_e_mva_nt_loose_2;
   Int_t           gen_match_2;
   Float_t         againstElectronLooseMVA6_2;
   Float_t         againstElectronMediumMVA6_2;
   Float_t         againstElectronTightMVA6_2;
   Float_t         againstElectronVLooseMVA6_2;
   Float_t         againstElectronVTightMVA6_2;
   Float_t         againstMuonLoose3_2;
   Float_t         againstMuonTight3_2;
   Float_t         byCombinedIsolationDeltaBetaCorrRaw3Hits_2;
   Float_t         byIsolationMVA3newDMwoLTraw_2;
   Float_t         byIsolationMVA3oldDMwoLTraw_2;
   Float_t         byIsolationMVA3newDMwLTraw_2;
   Float_t         byIsolationMVA3oldDMwLTraw_2;
   Float_t         chargedIsoPtSum_2;
   Float_t         decayModeFindingOldDMs_2;
   Float_t         neutralIsoPtSum_2;
   Float_t         puCorrPtSum_2;
   Float_t         trigweight_2;
   Float_t         idisoweight_2;
   Float_t         pt_tt;
   Float_t         mt_tot;
   Float_t         m_vis;
   Float_t         m_sv;
   Float_t         mt_sv;
   Float_t         met;
   Float_t         metphi;
   Float_t         puppimet;
   Float_t         puppimetphi;
   Float_t         mvamet;
   Float_t         mvametphi;
   Float_t         pzetavis;
   Float_t         pzetamiss;
   Float_t         pfpzetamiss;
   Float_t         puppipzetamiss;
   Float_t         mvacov00;
   Float_t         mvacov01;
   Float_t         mvacov10;
   Float_t         mvacov11;
   Float_t         metcov00;
   Float_t         metcov01;
   Float_t         metcov10;
   Float_t         metcov11;
   Float_t         mjj;
   Float_t         jdeta;
   Float_t         njetingap;
   Float_t         njetingap20;
   Float_t         jdphi;
   Float_t         nbtag;
   Float_t         njets;
   Float_t         njetspt20;
   Float_t         jpt_1;
   Float_t         jeta_1;
   Float_t         jphi_1;
   Float_t         jrawf_1;
   Float_t         jmva_1;
   Float_t         jpt_2;
   Float_t         jeta_2;
   Float_t         jphi_2;
   Float_t         jrawf_2;
   Float_t         jmva_2;
   Float_t         bpt_1;
   Float_t         beta_1;
   Float_t         bphi_1;
   Float_t         brawf_1;
   Float_t         bmva_1;
   Float_t         bcsv_1;
   Float_t         bpt_2;
   Float_t         beta_2;
   Float_t         bphi_2;
   Float_t         brawf_2;
   Float_t         bmva_2;
   Float_t         bcsv_2;
   Bool_t          dilepton_veto;
   Bool_t          extraelec_veto;
   Bool_t          extramuon_veto;
   Float_t         puweight;

   // List of branches
   TBranch        *b_eventWeight;   //!
   TBranch        *b_run;   //!
   TBranch        *b_lumi;   //!
   TBranch        *b_evt;   //!
   TBranch        *b_npv;   //!
   TBranch        *b_npu;   //!
   TBranch        *b_rho;   //!
   TBranch        *b_pt_1;   //!
   TBranch        *b_phi_1;   //!
   TBranch        *b_eta_1;   //!
   TBranch        *b_m_1;   //!
   TBranch        *b_q_1;   //!
   TBranch        *b_d0_1;   //!
   TBranch        *b_dZ_1;   //!
   TBranch        *b_mt_1;   //!
   TBranch        *b_pfmt_1;   //!
   TBranch        *b_puppimt_1;   //!
   TBranch        *b_iso_1;   //!
   TBranch        *b_id_e_mva_nt_loose_1;   //!
   TBranch        *b_gen_match_1;   //!
   TBranch        *b_againstElectronLooseMVA6_1;   //!
   TBranch        *b_againstElectronMediumMVA6_1;   //!
   TBranch        *b_againstElectronTightMVA6_1;   //!
   TBranch        *b_againstElectronVLooseMVA6_1;   //!
   TBranch        *b_againstElectronVTightMVA6_1;   //!
   TBranch        *b_againstMuonLoose3_1;   //!
   TBranch        *b_againstMuonTight3_1;   //!
   TBranch        *b_byCombinedIsolationDeltaBetaCorrRaw3Hits_1;   //!
   TBranch        *b_byIsolationMVA3newDMwoLTraw_1;   //!
   TBranch        *b_byIsolationMVA3oldDMwoLTraw_1;   //!
   TBranch        *b_byIsolationMVA3newDMwLTraw_1;   //!
   TBranch        *b_byIsolationMVA3oldDMwLTraw_1;   //!
   TBranch        *b_chargedIsoPtSum_1;   //!
   TBranch        *b_decayModeFindingOldDMs_1;   //!
   TBranch        *b_neutralIsoPtSum_1;   //!
   TBranch        *b_puCorrPtSum_1;   //!
   TBranch        *b_trigweight_1;   //!
   TBranch        *b_idisoweight_1;   //!
   TBranch        *b_pt_2;   //!
   TBranch        *b_phi_2;   //!
   TBranch        *b_eta_2;   //!
   TBranch        *b_m_2;   //!
   TBranch        *b_q_2;   //!
   TBranch        *b_d0_2;   //!
   TBranch        *b_dZ_2;   //!
   TBranch        *b_mt_2;   //!
   TBranch        *b_pfmt_2;   //!
   TBranch        *b_puppimt_2;   //!
   TBranch        *b_iso_2;   //!
   TBranch        *b_id_e_mva_nt_loose_2;   //!
   TBranch        *b_gen_match_2;   //!
   TBranch        *b_againstElectronLooseMVA6_2;   //!
   TBranch        *b_againstElectronMediumMVA6_2;   //!
   TBranch        *b_againstElectronTightMVA6_2;   //!
   TBranch        *b_againstElectronVLooseMVA6_2;   //!
   TBranch        *b_againstElectronVTightMVA6_2;   //!
   TBranch        *b_againstMuonLoose3_2;   //!
   TBranch        *b_againstMuonTight3_2;   //!
   TBranch        *b_byCombinedIsolationDeltaBetaCorrRaw3Hits_2;   //!
   TBranch        *b_byIsolationMVA3newDMwoLTraw_2;   //!
   TBranch        *b_byIsolationMVA3oldDMwoLTraw_2;   //!
   TBranch        *b_byIsolationMVA3newDMwLTraw_2;   //!
   TBranch        *b_byIsolationMVA3oldDMwLTraw_2;   //!
   TBranch        *b_chargedIsoPtSum_2;   //!
   TBranch        *b_decayModeFindingOldDMs_2;   //!
   TBranch        *b_neutralIsoPtSum_2;   //!
   TBranch        *b_puCorrPtSum_2;   //!
   TBranch        *b_trigweight_2;   //!
   TBranch        *b_idisoweight_2;   //!
   TBranch        *b_pt_tt;   //!
   TBranch        *b_mt_tot;   //!
   TBranch        *b_m_vis;   //!
   TBranch        *b_m_sv;   //!
   TBranch        *b_mt_sv;   //!
   TBranch        *b_met;   //!
   TBranch        *b_metphi;   //!
   TBranch        *b_puppimet;   //!
   TBranch        *b_puppimetphi;   //!
   TBranch        *b_mvamet;   //!
   TBranch        *b_mvametphi;   //!
   TBranch        *b_pzetavis;   //!
   TBranch        *b_pzetamiss;   //!
   TBranch        *b_pfpzetamiss;   //!
   TBranch        *b_puppipzetamiss;   //!
   TBranch        *b_mvacov00;   //!
   TBranch        *b_mvacov01;   //!
   TBranch        *b_mvacov10;   //!
   TBranch        *b_mvacov11;   //!
   TBranch        *b_metcov00;   //!
   TBranch        *b_metcov01;   //!
   TBranch        *b_metcov10;   //!
   TBranch        *b_metcov11;   //!
   TBranch        *b_mjj;   //!
   TBranch        *b_jdeta;   //!
   TBranch        *b_njetingap;   //!
   TBranch        *b_njetingap20;   //!
   TBranch        *b_jdphi;   //!
   TBranch        *b_nbtag;   //!
   TBranch        *b_njets;   //!
   TBranch        *b_njetspt20;   //!
   TBranch        *b_jpt_1;   //!
   TBranch        *b_jeta_1;   //!
   TBranch        *b_jphi_1;   //!
   TBranch        *b_jrawf_1;   //!
   TBranch        *b_jmva_1;   //!
   TBranch        *b_jpt_2;   //!
   TBranch        *b_jeta_2;   //!
   TBranch        *b_jphi_2;   //!
   TBranch        *b_jrawf_2;   //!
   TBranch        *b_jmva_2;   //!
   TBranch        *b_bpt_1;   //!
   TBranch        *b_beta_1;   //!
   TBranch        *b_bphi_1;   //!
   TBranch        *b_brawf_1;   //!
   TBranch        *b_bmva_1;   //!
   TBranch        *b_bcsv_1;   //!
   TBranch        *b_bpt_2;   //!
   TBranch        *b_beta_2;   //!
   TBranch        *b_bphi_2;   //!
   TBranch        *b_brawf_2;   //!
   TBranch        *b_bmva_2;   //!
   TBranch        *b_bcsv_2;   //!
   TBranch        *b_dilepton_veto;   //!
   TBranch        *b_extraelec_veto;   //!
   TBranch        *b_extramuon_veto;   //!
   TBranch        *b_puweight;   //!

   Synch(TTree *tree=0, std::string suffix="");
   virtual ~Synch();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef Synch_cxx
Synch::Synch(TTree *tree, std::string suffix) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   suffix_ = suffix;
   if (tree == 0) {
      /*
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("RootAnalysis_SynchNTuple.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("RootAnalysis_SynchNTuple.root");
      }
      TDirectory * dir = (TDirectory*)f->Get("RootAnalysis_SynchNTuple.root:/Summary");
      dir->GetObject("tree",tree);
      */
      TFile *f1 = new TFile("/afs/cern.ch/user/a/apyskir/scratch/HTauTau/RootAnalysis/HTTAnalysis/RootAnalysis_SynchNTuple.root");
      TChain *inChain = (TChain*)f1->Get("Summary/tree");
      inChain->Add("/afs/cern.ch/user/m/mspanrin/public/sync2016/SYNCFILE_SUSYGluGluToHToTauTau_M-160_mt_spring16_reHLT_v4.root/TauCheck/");
      tree = inChain;

   }
   Init(tree);
}

Synch::~Synch()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t Synch::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t Synch::LoadTree(Long64_t entry)
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

void Synch::Init(TTree *tree)
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

   fChain->SetBranchAddress("eventWeight", &eventWeight, &b_eventWeight);
   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("lumi", &lumi, &b_lumi);
   fChain->SetBranchAddress("evt", &evt, &b_evt);
   fChain->SetBranchAddress("npv", &npv, &b_npv);
   fChain->SetBranchAddress("npu", &npu, &b_npu);
   fChain->SetBranchAddress("rho", &rho, &b_rho);
   fChain->SetBranchAddress("pt_1", &pt_1, &b_pt_1);
   fChain->SetBranchAddress("phi_1", &phi_1, &b_phi_1);
   fChain->SetBranchAddress("eta_1", &eta_1, &b_eta_1);
   fChain->SetBranchAddress("m_1", &m_1, &b_m_1);
   fChain->SetBranchAddress("q_1", &q_1, &b_q_1);
   fChain->SetBranchAddress("d0_1", &d0_1, &b_d0_1);
   fChain->SetBranchAddress("dZ_1", &dZ_1, &b_dZ_1);
   fChain->SetBranchAddress("mt_1", &mt_1, &b_mt_1);
   fChain->SetBranchAddress("pfmt_1", &pfmt_1, &b_pfmt_1);
   fChain->SetBranchAddress("puppimt_1", &puppimt_1, &b_puppimt_1);
   fChain->SetBranchAddress("iso_1", &iso_1, &b_iso_1);
   fChain->SetBranchAddress("id_e_mva_nt_loose_1", &id_e_mva_nt_loose_1, &b_id_e_mva_nt_loose_1);
   fChain->SetBranchAddress("gen_match_1", &gen_match_1, &b_gen_match_1);
   fChain->SetBranchAddress("againstElectronLooseMVA6_1", &againstElectronLooseMVA6_1, &b_againstElectronLooseMVA6_1);
   fChain->SetBranchAddress("againstElectronMediumMVA6_1", &againstElectronMediumMVA6_1, &b_againstElectronMediumMVA6_1);
   fChain->SetBranchAddress("againstElectronTightMVA6_1", &againstElectronTightMVA6_1, &b_againstElectronTightMVA6_1);
   fChain->SetBranchAddress("againstElectronVLooseMVA6_1", &againstElectronVLooseMVA6_1, &b_againstElectronVLooseMVA6_1);
   fChain->SetBranchAddress("againstElectronVTightMVA6_1", &againstElectronVTightMVA6_1, &b_againstElectronVTightMVA6_1);
   fChain->SetBranchAddress("againstMuonLoose3_1", &againstMuonLoose3_1, &b_againstMuonLoose3_1);
   fChain->SetBranchAddress("againstMuonTight3_1", &againstMuonTight3_1, &b_againstMuonTight3_1);
   fChain->SetBranchAddress("byCombinedIsolationDeltaBetaCorrRaw3Hits_1", &byCombinedIsolationDeltaBetaCorrRaw3Hits_1, &b_byCombinedIsolationDeltaBetaCorrRaw3Hits_1);
   fChain->SetBranchAddress("byIsolationMVA3newDMwoLTraw_1", &byIsolationMVA3newDMwoLTraw_1, &b_byIsolationMVA3newDMwoLTraw_1);
   fChain->SetBranchAddress("byIsolationMVA3oldDMwoLTraw_1", &byIsolationMVA3oldDMwoLTraw_1, &b_byIsolationMVA3oldDMwoLTraw_1);
   fChain->SetBranchAddress("byIsolationMVA3newDMwLTraw_1", &byIsolationMVA3newDMwLTraw_1, &b_byIsolationMVA3newDMwLTraw_1);
   fChain->SetBranchAddress("byIsolationMVA3oldDMwLTraw_1", &byIsolationMVA3oldDMwLTraw_1, &b_byIsolationMVA3oldDMwLTraw_1);
   fChain->SetBranchAddress("chargedIsoPtSum_1", &chargedIsoPtSum_1, &b_chargedIsoPtSum_1);
   fChain->SetBranchAddress("decayModeFindingOldDMs_1", &decayModeFindingOldDMs_1, &b_decayModeFindingOldDMs_1);
   fChain->SetBranchAddress("neutralIsoPtSum_1", &neutralIsoPtSum_1, &b_neutralIsoPtSum_1);
   fChain->SetBranchAddress("puCorrPtSum_1", &puCorrPtSum_1, &b_puCorrPtSum_1);
   fChain->SetBranchAddress("trigweight_1", &trigweight_1, &b_trigweight_1);
   fChain->SetBranchAddress("idisoweight_1", &idisoweight_1, &b_idisoweight_1);
   fChain->SetBranchAddress("pt_2", &pt_2, &b_pt_2);
   fChain->SetBranchAddress("phi_2", &phi_2, &b_phi_2);
   fChain->SetBranchAddress("eta_2", &eta_2, &b_eta_2);
   fChain->SetBranchAddress("m_2", &m_2, &b_m_2);
   fChain->SetBranchAddress("q_2", &q_2, &b_q_2);
   fChain->SetBranchAddress("d0_2", &d0_2, &b_d0_2);
   fChain->SetBranchAddress("dZ_2", &dZ_2, &b_dZ_2);
   fChain->SetBranchAddress("mt_2", &mt_2, &b_mt_2);
   fChain->SetBranchAddress("pfmt_2", &pfmt_2, &b_pfmt_2);
   fChain->SetBranchAddress("puppimt_2", &puppimt_2, &b_puppimt_2);
   fChain->SetBranchAddress("iso_2", &iso_2, &b_iso_2);
   fChain->SetBranchAddress("id_e_mva_nt_loose_2", &id_e_mva_nt_loose_2, &b_id_e_mva_nt_loose_2);
   fChain->SetBranchAddress("gen_match_2", &gen_match_2, &b_gen_match_2);
   fChain->SetBranchAddress("againstElectronLooseMVA6_2", &againstElectronLooseMVA6_2, &b_againstElectronLooseMVA6_2);
   fChain->SetBranchAddress("againstElectronMediumMVA6_2", &againstElectronMediumMVA6_2, &b_againstElectronMediumMVA6_2);
   fChain->SetBranchAddress("againstElectronTightMVA6_2", &againstElectronTightMVA6_2, &b_againstElectronTightMVA6_2);
   fChain->SetBranchAddress("againstElectronVLooseMVA6_2", &againstElectronVLooseMVA6_2, &b_againstElectronVLooseMVA6_2);
   fChain->SetBranchAddress("againstElectronVTightMVA6_2", &againstElectronVTightMVA6_2, &b_againstElectronVTightMVA6_2);
   fChain->SetBranchAddress("againstMuonLoose3_2", &againstMuonLoose3_2, &b_againstMuonLoose3_2);
   fChain->SetBranchAddress("againstMuonTight3_2", &againstMuonTight3_2, &b_againstMuonTight3_2);
   fChain->SetBranchAddress("byCombinedIsolationDeltaBetaCorrRaw3Hits_2", &byCombinedIsolationDeltaBetaCorrRaw3Hits_2, &b_byCombinedIsolationDeltaBetaCorrRaw3Hits_2);
   fChain->SetBranchAddress("byIsolationMVA3newDMwoLTraw_2", &byIsolationMVA3newDMwoLTraw_2, &b_byIsolationMVA3newDMwoLTraw_2);
   fChain->SetBranchAddress("byIsolationMVA3oldDMwoLTraw_2", &byIsolationMVA3oldDMwoLTraw_2, &b_byIsolationMVA3oldDMwoLTraw_2);
   fChain->SetBranchAddress("byIsolationMVA3newDMwLTraw_2", &byIsolationMVA3newDMwLTraw_2, &b_byIsolationMVA3newDMwLTraw_2);
   fChain->SetBranchAddress("byIsolationMVA3oldDMwLTraw_2", &byIsolationMVA3oldDMwLTraw_2, &b_byIsolationMVA3oldDMwLTraw_2);
   fChain->SetBranchAddress("chargedIsoPtSum_2", &chargedIsoPtSum_2, &b_chargedIsoPtSum_2);
   fChain->SetBranchAddress("decayModeFindingOldDMs_2", &decayModeFindingOldDMs_2, &b_decayModeFindingOldDMs_2);
   fChain->SetBranchAddress("neutralIsoPtSum_2", &neutralIsoPtSum_2, &b_neutralIsoPtSum_2);
   fChain->SetBranchAddress("puCorrPtSum_2", &puCorrPtSum_2, &b_puCorrPtSum_2);
   fChain->SetBranchAddress("trigweight_2", &trigweight_2, &b_trigweight_2);
   fChain->SetBranchAddress("idisoweight_2", &idisoweight_2, &b_idisoweight_2);
   fChain->SetBranchAddress("pt_tt", &pt_tt, &b_pt_tt);
   fChain->SetBranchAddress("mt_tot", &mt_tot, &b_mt_tot);
   fChain->SetBranchAddress("m_vis", &m_vis, &b_m_vis);
   fChain->SetBranchAddress("m_sv", &m_sv, &b_m_sv);
   fChain->SetBranchAddress("mt_sv", &mt_sv, &b_mt_sv);
   fChain->SetBranchAddress("met", &met, &b_met);
   fChain->SetBranchAddress("metphi", &metphi, &b_metphi);
   fChain->SetBranchAddress("puppimet", &puppimet, &b_puppimet);
   fChain->SetBranchAddress("puppimetphi", &puppimetphi, &b_puppimetphi);
   fChain->SetBranchAddress("mvamet", &mvamet, &b_mvamet);
   fChain->SetBranchAddress("mvametphi", &mvametphi, &b_mvametphi);
   fChain->SetBranchAddress("pzetavis", &pzetavis, &b_pzetavis);
   fChain->SetBranchAddress("pzetamiss", &pzetamiss, &b_pzetamiss);
   fChain->SetBranchAddress("pfpzetamiss", &pfpzetamiss, &b_pfpzetamiss);
   fChain->SetBranchAddress("puppipzetamiss", &puppipzetamiss, &b_puppipzetamiss);
   fChain->SetBranchAddress("mvacov00", &mvacov00, &b_mvacov00);
   fChain->SetBranchAddress("mvacov01", &mvacov01, &b_mvacov01);
   fChain->SetBranchAddress("mvacov10", &mvacov10, &b_mvacov10);
   fChain->SetBranchAddress("mvacov11", &mvacov11, &b_mvacov11);
   fChain->SetBranchAddress("metcov00", &metcov00, &b_metcov00);
   fChain->SetBranchAddress("metcov01", &metcov01, &b_metcov01);
   fChain->SetBranchAddress("metcov10", &metcov10, &b_metcov10);
   fChain->SetBranchAddress("metcov11", &metcov11, &b_metcov11);
   fChain->SetBranchAddress("mjj", &mjj, &b_mjj);
   fChain->SetBranchAddress("jdeta", &jdeta, &b_jdeta);
   fChain->SetBranchAddress("njetingap", &njetingap, &b_njetingap);
   fChain->SetBranchAddress("njetingap20", &njetingap20, &b_njetingap20);
   fChain->SetBranchAddress("jdphi", &jdphi, &b_jdphi);
   fChain->SetBranchAddress("nbtag", &nbtag, &b_nbtag);
   fChain->SetBranchAddress("njets", &njets, &b_njets);
   fChain->SetBranchAddress("njetspt20", &njetspt20, &b_njetspt20);
   fChain->SetBranchAddress("jpt_1", &jpt_1, &b_jpt_1);
   fChain->SetBranchAddress("jeta_1", &jeta_1, &b_jeta_1);
   fChain->SetBranchAddress("jphi_1", &jphi_1, &b_jphi_1);
   fChain->SetBranchAddress("jrawf_1", &jrawf_1, &b_jrawf_1);
   fChain->SetBranchAddress("jmva_1", &jmva_1, &b_jmva_1);
   fChain->SetBranchAddress("jpt_2", &jpt_2, &b_jpt_2);
   fChain->SetBranchAddress("jeta_2", &jeta_2, &b_jeta_2);
   fChain->SetBranchAddress("jphi_2", &jphi_2, &b_jphi_2);
   fChain->SetBranchAddress("jrawf_2", &jrawf_2, &b_jrawf_2);
   fChain->SetBranchAddress("jmva_2", &jmva_2, &b_jmva_2);
   fChain->SetBranchAddress("bpt_1", &bpt_1, &b_bpt_1);
   fChain->SetBranchAddress("beta_1", &beta_1, &b_beta_1);
   fChain->SetBranchAddress("bphi_1", &bphi_1, &b_bphi_1);
   fChain->SetBranchAddress("brawf_1", &brawf_1, &b_brawf_1);
   fChain->SetBranchAddress("bmva_1", &bmva_1, &b_bmva_1);
   fChain->SetBranchAddress("bcsv_1", &bcsv_1, &b_bcsv_1);
   fChain->SetBranchAddress("bpt_2", &bpt_2, &b_bpt_2);
   fChain->SetBranchAddress("beta_2", &beta_2, &b_beta_2);
   fChain->SetBranchAddress("bphi_2", &bphi_2, &b_bphi_2);
   fChain->SetBranchAddress("brawf_2", &brawf_2, &b_brawf_2);
   fChain->SetBranchAddress("bmva_2", &bmva_2, &b_bmva_2);
   fChain->SetBranchAddress("bcsv_2", &bcsv_2, &b_bcsv_2);
   fChain->SetBranchAddress("dilepton_veto", &dilepton_veto, &b_dilepton_veto);
   fChain->SetBranchAddress("extraelec_veto", &extraelec_veto, &b_extraelec_veto);
   fChain->SetBranchAddress("extramuon_veto", &extramuon_veto, &b_extramuon_veto);
   fChain->SetBranchAddress("puweight", &puweight, &b_puweight);
   
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

Bool_t Synch::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void Synch::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t Synch::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef Synch_cxx
