/*#include "fastjet/ClusterSequence.hh"
#include "fastjet/PseudoJet.hh"
#include "TH1F.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include <vector>
#include <cmath>
#include <fstream>
#include <iostream>

using namespace fastjet;
using namespace std;

// Define the Matrix2x2 structure
struct Matrix2x2 {
  double m[2][2];
  Matrix2x2() {
    m[0][0] = m[0][1] = m[1][0] = m[1][1] = 0.0;
  }
};

// Function to calculate SL matrix for a single event
Matrix2x2 calculateSL(const vector<double>& px, const vector<double>& py) {
  Matrix2x2 SL;
  double sum_pt = 0.0; // Sum of transverse momenta
  
  for (size_t i = 0; i < px.size(); ++i) {
    double pt = sqrt(px[i] * px[i] + py[i] * py[i]);
    sum_pt += pt;
    
    if (pt > 0) {
      SL.m[0][0] += px[i] * px[i] / pt;
      SL.m[0][1] += px[i] * py[i] / pt;
      SL.m[1][0] += py[i] * px[i] / pt;
      SL.m[1][1] += py[i] * py[i] / pt;
    }
  }
  // Normalize SL by the sum of transverse momenta
  
  SL.m[0][0] /= sum_pt;
  SL.m[0][1] /= sum_pt;
  SL.m[1][0] /= sum_pt;
  SL.m[1][1] /= sum_pt;
  return SL;
  
}

// Function to calculate the eigenvalues of a 2x2 matrix
void calculateEigenvalues(const Matrix2x2& matrix, double& lambda1, double& lambda2) {
  double trace = matrix.m[0][0] + matrix.m[1][1];
  double det = matrix.m[0][0] * matrix.m[1][1] - matrix.m[0][1] * matrix.m[1][0];
  double discriminant = trace * trace - 4 * det;
  if (discriminant >= 0) {
    lambda1 = (trace + sqrt(discriminant)) / 2;
    lambda2 = (trace - sqrt(discriminant)) / 2;
  }
}

// Function to calculate transverse sphericity S_T
double calculateTransverseSphericity(const vector<double>& px, const vector<double>& py) {
  Matrix2x2 SL = calculateSL(px, py);
  double lambda1, lambda2;
  calculateEigenvalues(SL, lambda1, lambda2);
  
  // Sort lambda1 and lambda2 to ensure lambda1 is always greater or equal
  if (lambda1 < lambda2) swap(lambda1, lambda2);
  
  return 2 * lambda2 / (lambda1 + lambda2);
}



void subjets() {
  // Open the ROOT file and get the tree
  TFile *inputFile = TFile::Open("/home/shalini/Software/pythia/examples/13mpion_shalini.root");
  TTree *tree = (TTree *)inputFile->Get("Tree");
  
  // Set branch addresses to read particle data
  Int_t numParticles;
  Double_t px[500], py[500], pz[500], energy[500], pt[500], eta[500], phi[500];
  tree->SetBranchAddress("numParticles", &numParticles);
  tree->SetBranchAddress("px", px);
  tree->SetBranchAddress("py", py);
  tree->SetBranchAddress("pz", pz);
  tree->SetBranchAddress("energy", energy);
  tree->SetBranchAddress("pt", pt);
  tree->SetBranchAddress("eta", eta);
  tree->SetBranchAddress("phi", phi);
  
  // Histograms for transverse sphericity and number of jets
  TH1F *sphericityHist_AllJets = new TH1F("sphericityHist_AllJets", "Transverse Sphericity;S_{T};dN/dS_{T}", 100, 0, 1);
  TH1F *sphericityHist_0Jets = new TH1F("sphericityHist_0Jets", "Transverse Sphericity (0 Jets);S_{T};dN/dS_{T}", 100, 0, 1);
  TH1F *sphericityHist_1Jet = new TH1F("sphericityHist_1Jet", "Transverse Sphericity (1 Jet);S_{T};dN/dS_{T}", 100, 0, 1);
  TH1F *sphericityHist_2Jets = new TH1F("sphericityHist_2Jets", "Transverse Sphericity (2 Jets);S_{T};dN/dS_{T}", 100, 0, 1);
  TH1F *sphericityHist_3PlusJets = new TH1F("sphericityHist_3PlusJets", "Transverse Sphericity (3+ Jets);S_{T};dN/dS_{T}", 100, 0, 1);
  
  // Jet clustering parameters
  double EtMin = 10;
  double R = 0.4;
  JetDefinition jet_def(antikt_algorithm, R);
  vector<PseudoJet> particles;
  
  Int_t nentries = (Int_t)tree->GetEntries();
  
  // Open the log file
  ofstream logFile("debug_log.txt");
  
  // Counter for valid entries
  int validEntries = 0;
  
  // Add a counter for events with 0 jets and S_T = 0
  int zeroJets_ST_zero_count = 0;
  
  for (Int_t i = 0; i < nentries; i++) {
    tree->GetEntry(i);
    particles.clear();
    vector<double> event_px, event_py;
    
    // Prepare particles for jet clustering and sphericity calculation
    for (Int_t j = 0; j < numParticles; j++) {
      particles.push_back(PseudoJet(px[j], py[j], pz[j], energy[j]));
      event_px.push_back(px[j]);
      event_py.push_back(py[j]);
    }
    
    // Cluster jets
    ClusterSequence clustSeq(particles, jet_def);
    Selector jet_selector = SelectorEtMin(EtMin);
    vector<PseudoJet> inclusive_jets = sorted_by_pt(jet_selector(clustSeq.inclusive_jets()));
    
    // Count the number of jets
    int numJets = inclusive_jets.size();
    
    // Calculate transverse sphericity S_T for the event
    double ST = calculateTransverseSphericity(event_px, event_py);
    
    // Log the number of jets and ST value
    logFile << "Event " << i << ": numJets = " << numJets << ", S_T = " << ST << endl;
    
    // Only fill the histograms if S_T is valid (between 0 and 1)
    if (ST >= 0 && ST <= 1) {
      validEntries++;
      // Fill the AllJets histogram
      sphericityHist_AllJets->Fill(ST);
      
      // Fill the corresponding jet multiplicity histogram
      if (numJets < 1) {
	sphericityHist_0Jets->Fill(ST);
	if (ST == 0) {
	  zeroJets_ST_zero_count++;
	}
      } else if (numJets == 1) {
	sphericityHist_1Jet->Fill(ST);
      } else if (numJets == 2) {
	sphericityHist_2Jets->Fill(ST);
      } else if (numJets >= 3) {
	sphericityHist_3PlusJets->Fill(ST);
      }
    } else {
      logFile << "Event " << i << " has invalid S_T value: " << ST << endl;
    }
  }
  
  // Print the number of events with 0 jets and S_T = 0
  cout << "Number of events with 0 jets and S_T = 0: " << zeroJets_ST_zero_count << endl;
  
  
  // Close the log file
  logFile.close();
  
  // Normalize histograms
  double scale = 1.0 / validEntries;
  sphericityHist_AllJets->Scale(scale);
  sphericityHist_0Jets->Scale(scale);
  sphericityHist_1Jet->Scale(scale);
  sphericityHist_2Jets->Scale(scale);
  sphericityHist_3PlusJets->Scale(scale);
  
  // Function to style the histogram and canvas
  auto styleHistogram = [](TH1F* hist, Color_t color) {
    hist->SetMarkerStyle(21);
    hist->SetMarkerColor(color);
    hist->SetLineColor(color);
    hist->GetXaxis()->SetTitleSize(0.04);
    hist->GetYaxis()->SetTitleSize(0.04);
    hist->SetStats(1); // Display the statistics box
  };
  
  // Style and draw histograms for each jet multiplicity
  styleHistogram(sphericityHist_AllJets, kMagenta);
  styleHistogram(sphericityHist_0Jets, kBlack);
  styleHistogram(sphericityHist_1Jet, kRed);
  styleHistogram(sphericityHist_2Jets, kBlue);
  styleHistogram(sphericityHist_3PlusJets, kGreen);
  
  
  
  TCanvas *c1 = new TCanvas("c1", "Transverse Sphericity - All Jets", 800, 600);
  c1->SetLogy(); // Log scale on y-axis
  sphericityHist_AllJets->Draw("P E"); // Draw points with error bars
  //c1->SaveAs("13.6mpion_TransverseSphericity_AllJets.png");
  
  TCanvas *c2 = new TCanvas("c2", "Transverse Sphericity - 0 Jets", 800, 600);
  c2->SetLogy();
  sphericityHist_0Jets->Draw("P E");
  //c2->SaveAs("13.6mpion_TransverseSphericity_0Jets.png");
  
  TCanvas *c3 = new TCanvas("c3", "Transverse Sphericity - 1 Jet", 800, 600);
  c3->SetLogy();
  sphericityHist_1Jet->Draw("P E");
  //c3->SaveAs("13.6mpion_TransverseSphericity_1Jet.png");

  
  TCanvas *c4 = new TCanvas("c4", "Transverse Sphericity - 2 Jets", 800, 600);
  c4->SetLogy();
  sphericityHist_2Jets->Draw("P E");
  //c4->SaveAs("13.6mpion_TransverseSphericity_2Jets.png");

  
  TCanvas *c5 = new TCanvas("c5", "Transverse Sphericity - 3+ Jets", 800, 600);
  c5->SetLogy();
  sphericityHist_3PlusJets->Draw("P E");
  //c5->SaveAs("13.6mpion_TransverseSphericity_3PlusJets.png");
  
  // Save histograms to a single ROOT file
  TFile *outputFile = new TFile("jetplots_13mpion.root", "RECREATE");
  sphericityHist_AllJets->Write();
  sphericityHist_0Jets->Write();
  sphericityHist_1Jet->Write();
  sphericityHist_2Jets->Write();
  sphericityHist_3PlusJets->Write();
  //outputFile->Close();
 
  
}


*/

#include "fastjet/ClusterSequence.hh"
#include "fastjet/PseudoJet.hh"
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TPaveStats.h"
#include <vector>
#include <cmath>
#include <fstream>
#include <iostream>
#include <iomanip>

using namespace fastjet;
using namespace std;

// [Previous Matrix2x2 struct and helper functions remain exactly the same]
struct Matrix2x2 {
    double m[2][2];
    Matrix2x2() {
        m[0][0] = m[0][1] = m[1][0] = m[1][1] = 0.0;
    }
};

Matrix2x2 calculateSL(const vector<double>& px, const vector<double>& py) {
    Matrix2x2 SL;
    double sum_pt = 0.0;
    
    for (size_t i = 0; i < px.size(); ++i) {
        double pt = sqrt(px[i] * px[i] + py[i] * py[i]);
        sum_pt += pt;
        
        if (pt > 0) {
            SL.m[0][0] += px[i] * px[i] / pt;
            SL.m[0][1] += px[i] * py[i] / pt;
            SL.m[1][0] += py[i] * px[i] / pt;
            SL.m[1][1] += py[i] * py[i] / pt;
        }
    }
    
    if (sum_pt > 0) {
        SL.m[0][0] /= sum_pt;
        SL.m[0][1] /= sum_pt;
        SL.m[1][0] /= sum_pt;
        SL.m[1][1] /= sum_pt;
    }
    
    return SL;
}

void calculateEigenvalues(const Matrix2x2& matrix, double& lambda1, double& lambda2) {
    double trace = matrix.m[0][0] + matrix.m[1][1];
    double det = matrix.m[0][0] * matrix.m[1][1] - matrix.m[0][1] * matrix.m[1][0];
    double discriminant = trace * trace - 4 * det;
    
    if (discriminant >= 0) {
        lambda1 = (trace + sqrt(discriminant)) / 2;
        lambda2 = (trace - sqrt(discriminant)) / 2;
    } else {
        lambda1 = lambda2 = trace / 2;
    }
}

double calculateTransverseSphericity(const vector<double>& px, const vector<double>& py) {
    if (px.empty() || py.empty()) {
        return 0.0;
    }
    
    Matrix2x2 SL = calculateSL(px, py);
    double lambda1, lambda2;
    calculateEigenvalues(SL, lambda1, lambda2);
    
    if (lambda1 < lambda2) {
        swap(lambda1, lambda2);
    }
    
    if ((lambda1 + lambda2) > 0) {
        return 2 * lambda2 / (lambda1 + lambda2);
    }
    return 0.0;
}

void subjets() {
    // Open the ROOT file and get the tree
    TFile *inputFile = TFile::Open("/home/shalini/Software/pythia/examples/13mpion_shalini.root");
    if (!inputFile || inputFile->IsZombie()) {
        cout << "Error opening input file!" << endl;
        return;
    }
    
    TTree *tree = (TTree *)inputFile->Get("Tree");
    if (!tree) {
        cout << "Error accessing tree!" << endl;
        inputFile->Close();
        return;
    }
    
    // Set branch addresses
    Int_t numParticles;
    Double_t px[500], py[500], pz[500], energy[500], pt[500], eta[500], phi[500];
    tree->SetBranchAddress("numParticles", &numParticles);
    tree->SetBranchAddress("px", px);
    tree->SetBranchAddress("py", py);
    tree->SetBranchAddress("pz", pz);
    tree->SetBranchAddress("energy", energy);
    tree->SetBranchAddress("pt", pt);
    tree->SetBranchAddress("eta", eta);
    tree->SetBranchAddress("phi", phi);
    
    // Create histograms
    const int nBins = 50;
    const double xMin = 0.0;
    const double xMax = 1.0;
    
    TH1F *sphericityHist_AllJets = new TH1F("sphericityHist_AllJets", "Transverse Sphericity;S_{T};(1/N_{ev})dN/dS_{T}", nBins, xMin, xMax);
    TH1F *sphericityHist_0Jets = new TH1F("sphericityHist_0Jets", "Transverse Sphericity (0 Jets);S_{T};(1/N_{ev})dN/dS_{T}", nBins, xMin, xMax);
    TH1F *sphericityHist_1Jet = new TH1F("sphericityHist_1Jet", "Transverse Sphericity (1 Jet);S_{T};(1/N_{ev})dN/dS_{T}", nBins, xMin, xMax);
    TH1F *sphericityHist_2Jets = new TH1F("sphericityHist_2Jets", "Transverse Sphericity (2 Jets);S_{T};(1/N_{ev})dN/dS_{T}", nBins, xMin, xMax);
    TH1F *sphericityHist_3PlusJets = new TH1F("sphericityHist_3PlusJets", "Transverse Sphericity (3+ Jets);S_{T};(1/N_{ev})dN/dS_{T}", nBins, xMin, xMax);
    
    // Diagnostic histograms
    TH1F *rawJetMultiplicity = new TH1F("rawJetMultiplicity", "Jet Multiplicity (No Et Cut);N_{jets};Events", 20, -0.5, 19.5);
    TH1F *filteredJetMultiplicity = new TH1F("filteredJetMultiplicity", "Jet Multiplicity (Et > 10 GeV);N_{jets};Events", 20, -0.5, 19.5);
    TH1F *eventEtDist = new TH1F("eventEtDist", "Total Event Et;E_{T} [GeV];Events", 100, 0, 500);
    TH2F *jetEtVsST = new TH2F("jetEtVsST", "Leading Jet Et vs S_{T};S_{T};Leading Jet E_{T} [GeV]", nBins, xMin, xMax, 100, 0, 200);
    
    // Jet clustering parameters
    double EtMin = 10.0; // GeV
    double R = 0.4;
    JetDefinition jet_def(antikt_algorithm, R);
    vector<PseudoJet> particles;
    
    Int_t nentries = (Int_t)tree->GetEntries();
    
    // Event counting variables
    int validEntries = 0;
    int zeroJetEvents = 0;
    int oneJetEvents = 0;
    int twoJetEvents = 0;
    int threePlusJetEvents = 0;
    
    // Main event loop
    for (Int_t i = 0; i < nentries; i++) {
        if (i % 1000 == 0) {
            cout << "Processing event " << i << "/" << nentries << "\r" << flush;
        }
        
        tree->GetEntry(i);
        particles.clear();
        vector<double> event_px, event_py;
        double total_event_Et = 0.0;
        
        for (Int_t j = 0; j < numParticles; j++) {
            if (std::isfinite(px[j]) && std::isfinite(py[j]) && 
                std::isfinite(pz[j]) && std::isfinite(energy[j])) {
                particles.push_back(PseudoJet(px[j], py[j], pz[j], energy[j]));
                event_px.push_back(px[j]);
                event_py.push_back(py[j]);
                total_event_Et += sqrt(px[j]*px[j] + py[j]*py[j]);
            }
        }
        
        if (particles.empty()) continue;
        
        // Cluster jets
        ClusterSequence clustSeq(particles, jet_def);
        vector<PseudoJet> all_jets = sorted_by_pt(clustSeq.inclusive_jets());
        rawJetMultiplicity->Fill(all_jets.size());
        
        // Apply Et cut
        Selector jet_selector = SelectorEtMin(EtMin);
        vector<PseudoJet> inclusive_jets = sorted_by_pt(jet_selector(all_jets));
        filteredJetMultiplicity->Fill(inclusive_jets.size());
        
        // Calculate transverse sphericity
        double ST = calculateTransverseSphericity(event_px, event_py);
        
        if (ST >= 0 && ST <= 1) {
            validEntries++;
            double weight = 1.0;
            
            sphericityHist_AllJets->Fill(ST, weight);
            
            if (inclusive_jets.empty()) {
                sphericityHist_0Jets->Fill(ST, weight);
                zeroJetEvents++;
            } else if (inclusive_jets.size() == 1) {
                sphericityHist_1Jet->Fill(ST, weight);
                oneJetEvents++;
            } else if (inclusive_jets.size() == 2) {
                sphericityHist_2Jets->Fill(ST, weight);
                twoJetEvents++;
            } else {
                sphericityHist_3PlusJets->Fill(ST, weight);
                threePlusJetEvents++;
            }
            
            eventEtDist->Fill(total_event_Et);
            if (!all_jets.empty()) {
                jetEtVsST->Fill(ST, all_jets[0].Et());
            }
        }
    }
    
    // Normalize histograms
    if (validEntries > 0) {
        double binWidth = (xMax - xMin) / nBins;
        double normFactor = 1.0 / (validEntries * binWidth);
        
        sphericityHist_AllJets->Scale(normFactor);
        sphericityHist_0Jets->Scale(normFactor);
        sphericityHist_1Jet->Scale(normFactor);
        sphericityHist_2Jets->Scale(normFactor);
        sphericityHist_3PlusJets->Scale(normFactor);
    }
    
    // Find global maximum and minimum across all histograms
    double globalMax = -1e10;
    double globalMin = 1e10;
    
    TH1F* allHists[] = {
        sphericityHist_AllJets,
        sphericityHist_0Jets,
        sphericityHist_1Jet,
        sphericityHist_2Jets,
        sphericityHist_3PlusJets
    };
    
    for (auto hist : allHists) {
        for (int bin = 1; bin <= hist->GetNbinsX(); bin++) {
            double content = hist->GetBinContent(bin);
            if (content > 0) {  // Only consider positive values for log scale
                globalMax = max(globalMax, content);
                globalMin = min(globalMin, content);
            }
        }
    }
    
    // Add some padding to the ranges
    double logMin = log10(globalMin);
    double logMax = log10(globalMax);
    double logRange = logMax - logMin;
    
    double yMin = pow(10, logMin - 0.1 * logRange);
    double yMax = pow(10, logMax + 0.1 * logRange);
    
    // Style histograms
    auto styleHistogram = [](TH1F* hist, Color_t color) {
        hist->SetMarkerStyle(20);
        hist->SetMarkerSize(1.2);
        hist->SetMarkerColor(color);
        hist->SetLineColor(color);
        hist->SetLineWidth(2);
        hist->GetXaxis()->SetTitleSize(0.045);
        hist->GetYaxis()->SetTitleSize(0.045);
        hist->GetXaxis()->SetLabelSize(0.04);
        hist->GetYaxis()->SetLabelSize(0.04);
    };
    
    styleHistogram(sphericityHist_AllJets, kBlack);
    styleHistogram(sphericityHist_0Jets, kRed);
    styleHistogram(sphericityHist_1Jet, kBlue);
    styleHistogram(sphericityHist_2Jets, kGreen+2);
    styleHistogram(sphericityHist_3PlusJets, kMagenta);
    
    // Set y-axis range for all histograms
    for (auto hist : allHists) {
        hist->GetYaxis()->SetRangeUser(yMin, yMax);
    }
    
    // Create combined plot
    TCanvas *cCombined = new TCanvas("cCombined", "Combined Sphericity", 800, 600);
    cCombined->SetLogy();
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.05);
    
    sphericityHist_AllJets->Draw("P E");
    sphericityHist_0Jets->Draw("P E SAME");
    sphericityHist_1Jet->Draw("P E SAME");
    sphericityHist_2Jets->Draw("P E SAME");
    sphericityHist_3PlusJets->Draw("P E SAME");
    
    // Create legend
    TLegend *legend = new TLegend(0.65, 0.65, 0.85, 0.85);
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    legend->AddEntry(sphericityHist_AllJets, "All Jets", "p");
    legend->AddEntry(sphericityHist_0Jets, "0 Jets", "p");
    legend->AddEntry(sphericityHist_1Jet, "1 Jet", "p");
    legend->AddEntry(sphericityHist_2Jets, "2 Jets", "p");
    legend->AddEntry(sphericityHist_3PlusJets, "3+ Jets", "p");
    legend->Draw();

    // Add event statistics
    TPaveText *stats = new TPaveText(0.65, 0.45, 0.85, 0.63, "NDC");
    stats->SetBorderSize(0);
    stats->SetFillStyle(0);
    stats->AddText(Form("Total Events: %d", validEntries));
    stats->AddText(Form("0 Jets: %d", zeroJetEvents));
    stats->AddText(Form("1 Jet: %d", oneJetEvents));
    stats->AddText(Form("2 Jets: %d", twoJetEvents));
    stats->AddText(Form("3+ Jets: %d", threePlusJetEvents));
    stats->Draw();
    
    // Create individual canvases with consistent y-axis ranges
    TCanvas *c1 = new TCanvas("c1", "Transverse Sphericity - All Jets", 800, 600);
    c1->SetLogy();
    gPad->SetLeftMargin(0.15);
    sphericityHist_AllJets->Draw("P E");
    sphericityHist_AllJets->GetYaxis()->SetRangeUser(yMin, yMax);
    
    TCanvas *c2 = new TCanvas("c2", "Transverse Sphericity - 0 Jets", 800, 600);
    c2->SetLogy();
    gPad->SetLeftMargin(0.15);
    sphericityHist_0Jets->Draw("P E");
    sphericityHist_0Jets->GetYaxis()->SetRangeUser(yMin, yMax);
    
    TCanvas *c3 = new TCanvas("c3", "Transverse Sphericity - 1 Jet", 800, 600);
    c3->SetLogy();
    gPad->SetLeftMargin(0.15);
    sphericityHist_1Jet->Draw("P E");
    sphericityHist_1Jet->GetYaxis()->SetRangeUser(yMin, yMax);
    
    TCanvas *c4 = new TCanvas("c4", "Transverse Sphericity - 2 Jets", 800, 600);
    c4->SetLogy();
    gPad->SetLeftMargin(0.15);
    sphericityHist_2Jets->Draw("P E");
    sphericityHist_2Jets->GetYaxis()->SetRangeUser(yMin, yMax);
    
    TCanvas *c5 = new TCanvas("c5", "Transverse Sphericity - 3+ Jets", 800, 600);
    c5->SetLogy();
    gPad->SetLeftMargin(0.15);
    sphericityHist_3PlusJets->Draw("P E");
    sphericityHist_3PlusJets->GetYaxis()->SetRangeUser(yMin, yMax);
    
    // Diagnostic plots
    TCanvas *cJetMult = new TCanvas("cJetMult", "Jet Multiplicity Comparison", 800, 800);
    cJetMult->Divide(1, 2);
    
    cJetMult->cd(1);
    gPad->SetLogy();
    rawJetMultiplicity->Draw();
    
    cJetMult->cd(2);
    gPad->SetLogy();
    filteredJetMultiplicity->Draw();
    
    TCanvas *cEventEt = new TCanvas("cEventEt", "Event Et Distribution", 800, 600);
    cEventEt->SetLogy();
    eventEtDist->Draw();
    
    TCanvas *cJetEtVsST = new TCanvas("cJetEtVsST", "Leading Jet Et vs ST", 800, 600);
    jetEtVsST->Draw("COLZ");
    
    // Save all plots to ROOT file
    TFile *outputFile = new TFile("jetplots_13mpion_improved.root", "RECREATE");
    
    // Write histograms
    sphericityHist_AllJets->Write();
    sphericityHist_0Jets->Write();
    sphericityHist_1Jet->Write();
    sphericityHist_2Jets->Write();
    sphericityHist_3PlusJets->Write();
    rawJetMultiplicity->Write();
    filteredJetMultiplicity->Write();
    eventEtDist->Write();
    jetEtVsST->Write();
    
    // Write canvases
    cCombined->Write();
    c1->Write();
    c2->Write();
    c3->Write();
    c4->Write();
    c5->Write();
    cJetMult->Write();
    cEventEt->Write();
    cJetEtVsST->Write();
    
    // Save as images
    //cCombined->SaveAs("13mpion_TransverseSphericity_Combined.png");
    //c1->SaveAs("13mpion_TransverseSphericity_AllJets.png");
    //c2->SaveAs("13mpion_TransverseSphericity_0Jets.png");
    //c3->SaveAs("13mpion_TransverseSphericity_1Jet.png");
    //c4->SaveAs("13mpion_TransverseSphericity_2Jets.png");
    //c5->SaveAs("13mpion_TransverseSphericity_3PlusJets.png");
    //cJetMult->SaveAs("13mpion_JetMultiplicity.png");
    //cEventEt->SaveAs("13mpion_EventEt.png");
    //cJetEtVsST->SaveAs("13mpion_JetEtVsST.png");
    
    // Clean up
    // outputFile->Close();
    // inputFile->Close();
    
    // Print summary
    //cout << "\nAnalysis completed successfully!" << endl;
    //cout << "Total events processed: " << nentries << endl;
    //cout << "Valid events: " << validEntries << endl;
    //cout << "Events by jet multiplicity:" << endl;
    //cout << "  0 jets: " << zeroJetEvents << endl;
    //cout << "  1 jet:  " << oneJetEvents << endl;
    //cout << "  2 jets: " << twoJetEvents << endl;
    //cout << "  3+ jets: " << threePlusJetEvents << endl;
    //cout << "Output saved to jetplots_13mpion_improved.root" << endl;
}







/*
#include "fastjet/ClusterSequence.hh"
#include "fastjet/PseudoJet.hh"
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TPaveStats.h"
#include <vector>
#include <cmath>
#include <fstream>
#include <iostream>
#include <iomanip>

using namespace fastjet;
using namespace std;

struct Matrix2x2 {
    double m[2][2];
    Matrix2x2() {
        m[0][0] = m[0][1] = m[1][0] = m[1][1] = 0.0;
    }
};

Matrix2x2 calculateSL(const vector<double>& px, const vector<double>& py) {
    Matrix2x2 SL;
    double sum_pt = 0.0;
    
    for (size_t i = 0; i < px.size(); ++i) {
        double pt = sqrt(px[i] * px[i] + py[i] * py[i]);
        sum_pt += pt;
        
        if (pt > 0) {
            SL.m[0][0] += px[i] * px[i] / pt;
            SL.m[0][1] += px[i] * py[i] / pt;
            SL.m[1][0] += py[i] * px[i] / pt;
            SL.m[1][1] += py[i] * py[i] / pt;
        }
    }
    
    if (sum_pt > 0) {
        SL.m[0][0] /= sum_pt;
        SL.m[0][1] /= sum_pt;
        SL.m[1][0] /= sum_pt;
        SL.m[1][1] /= sum_pt;
    }
    
    return SL;
}

void calculateEigenvalues(const Matrix2x2& matrix, double& lambda1, double& lambda2) {
    double trace = matrix.m[0][0] + matrix.m[1][1];
    double det = matrix.m[0][0] * matrix.m[1][1] - matrix.m[0][1] * matrix.m[1][0];
    double discriminant = trace * trace - 4 * det;
    
    if (discriminant >= 0) {
        lambda1 = (trace + sqrt(discriminant)) / 2;
        lambda2 = (trace - sqrt(discriminant)) / 2;
    } else {
        lambda1 = lambda2 = trace / 2;
    }
}

double calculateTransverseSphericity(const vector<double>& px, const vector<double>& py) {
    if (px.empty() || py.empty()) {
        return 0.0;
    }
    
    Matrix2x2 SL = calculateSL(px, py);
    double lambda1, lambda2;
    calculateEigenvalues(SL, lambda1, lambda2);
    
    if (lambda1 < lambda2) {
        swap(lambda1, lambda2);
    }
    
    if ((lambda1 + lambda2) > 0) {
        return 2 * lambda2 / (lambda1 + lambda2);
    }
    return 0.0;
}

void subjets() {
    // Open the ROOT file and get the tree
    TFile *inputFile = TFile::Open("/home/shalini/Software/pythia/examples/13mpion_shalini.root");
    if (!inputFile || inputFile->IsZombie()) {
        cout << "Error opening input file!" << endl;
        return;
    }
    
    TTree *tree = (TTree *)inputFile->Get("Tree");
    if (!tree) {
        cout << "Error accessing tree!" << endl;
        inputFile->Close();
        return;
    }
    
    // Set branch addresses
    Int_t numParticles;
    Double_t px[500], py[500], pz[500], energy[500], pt[500], eta[500], phi[500];
    tree->SetBranchAddress("numParticles", &numParticles);
    tree->SetBranchAddress("px", px);
    tree->SetBranchAddress("py", py);
    tree->SetBranchAddress("pz", pz);
    tree->SetBranchAddress("energy", energy);
    tree->SetBranchAddress("pt", pt);
    tree->SetBranchAddress("eta", eta);
    tree->SetBranchAddress("phi", phi);
    
    // Create histograms
    const int nBins = 50;
    const double xMin = 0.0;
    const double xMax = 1.0;
    
    TH1F *sphericityHist_AllJets = new TH1F("sphericityHist_AllJets", "Transverse Sphericity;S_{T};(1/N_{ev})dN/dS_{T}", nBins, xMin, xMax);
    TH1F *sphericityHist_0Jets = new TH1F("sphericityHist_0Jets", "Transverse Sphericity (0 Jets);S_{T};(1/N_{ev})dN/dS_{T}", nBins, xMin, xMax);
    TH1F *sphericityHist_1Jet = new TH1F("sphericityHist_1Jet", "Transverse Sphericity (1 Jet);S_{T};(1/N_{ev})dN/dS_{T}", nBins, xMin, xMax);
    TH1F *sphericityHist_2Jets = new TH1F("sphericityHist_2Jets", "Transverse Sphericity (2 Jets);S_{T};(1/N_{ev})dN/dS_{T}", nBins, xMin, xMax);
    TH1F *sphericityHist_3PlusJets = new TH1F("sphericityHist_3PlusJets", "Transverse Sphericity (3+ Jets);S_{T};(1/N_{ev})dN/dS_{T}", nBins, xMin, xMax);
    
    // Diagnostic histograms
    TH1F *rawJetMultiplicity = new TH1F("rawJetMultiplicity", "Jet Multiplicity (No Et Cut);N_{jets};Events", 20, -0.5, 19.5);
    TH1F *filteredJetMultiplicity = new TH1F("filteredJetMultiplicity", "Jet Multiplicity (Et > 10 GeV);N_{jets};Events", 20, -0.5, 19.5);
    TH1F *eventEtDist = new TH1F("eventEtDist", "Total Event Et;E_{T} [GeV];Events", 100, 0, 500);
    TH2F *jetEtVsST = new TH2F("jetEtVsST", "Leading Jet Et vs S_{T};S_{T};Leading Jet E_{T} [GeV]", nBins, xMin, xMax, 100, 0, 200);
    
    // Jet clustering parameters
    double EtMin = 10.0; // GeV
    double R = 0.4;
    JetDefinition jet_def(antikt_algorithm, R);
    vector<PseudoJet> particles;
    
    Int_t nentries = (Int_t)tree->GetEntries();
    
    // Event counting variables
    int validEntries = 0;
    int zeroJetEvents = 0;
    int oneJetEvents = 0;
    int twoJetEvents = 0;
    int threePlusJetEvents = 0;
    
    // Main event loop
    for (Int_t i = 0; i < nentries; i++) {
        if (i % 1000 == 0) {
            cout << "Processing event " << i << "/" << nentries << "\r" << flush;
        }
        
        tree->GetEntry(i);
        particles.clear();
        vector<double> event_px, event_py;
        double total_event_Et = 0.0;
        
        for (Int_t j = 0; j < numParticles; j++) {
            if (std::isfinite(px[j]) && std::isfinite(py[j]) && 
                std::isfinite(pz[j]) && std::isfinite(energy[j])) {
                particles.push_back(PseudoJet(px[j], py[j], pz[j], energy[j]));
                event_px.push_back(px[j]);
                event_py.push_back(py[j]);
                total_event_Et += sqrt(px[j]*px[j] + py[j]*py[j]);
            }
        }
        
        if (particles.empty()) continue;
        
        // Cluster jets
        ClusterSequence clustSeq(particles, jet_def);
        vector<PseudoJet> all_jets = sorted_by_pt(clustSeq.inclusive_jets());
        rawJetMultiplicity->Fill(all_jets.size());
        
        // Apply Et cut
        Selector jet_selector = SelectorEtMin(EtMin);
        vector<PseudoJet> inclusive_jets = sorted_by_pt(jet_selector(all_jets));
        filteredJetMultiplicity->Fill(inclusive_jets.size());
        
        // Calculate transverse sphericity
        double ST = calculateTransverseSphericity(event_px, event_py);
        
        if (ST >= 0 && ST <= 1) {
            validEntries++;
            double weight = 1.0;
            
            sphericityHist_AllJets->Fill(ST, weight);
            
            if (inclusive_jets.empty()) {
                sphericityHist_0Jets->Fill(ST, weight);
                zeroJetEvents++;
            } else if (inclusive_jets.size() == 1) {
                sphericityHist_1Jet->Fill(ST, weight);
                oneJetEvents++;
            } else if (inclusive_jets.size() == 2) {
                sphericityHist_2Jets->Fill(ST, weight);
                twoJetEvents++;
            } else {
                sphericityHist_3PlusJets->Fill(ST, weight);
                threePlusJetEvents++;
            }
            
            // Fill additional histograms
            eventEtDist->Fill(total_event_Et);
            if (!inclusive_jets.empty()) {
                jetEtVsST->Fill(ST, inclusive_jets[0].pt());
            }
        }
    }
    
    // Close input file
    inputFile->Close();
    
    // Normalize histograms
    if (validEntries > 0) {
        double binWidth = (xMax - xMin) / nBins;
        double normFactor = 1.0 / (validEntries * binWidth);
        
        sphericityHist_AllJets->Scale(normFactor);
        sphericityHist_0Jets->Scale(normFactor);
        sphericityHist_1Jet->Scale(normFactor);
        sphericityHist_2Jets->Scale(normFactor);
        sphericityHist_3PlusJets->Scale(normFactor);
    }

    // Find global maximum and minimum across all histograms
    double globalMax = -1e10;
    double globalMin = 1e10;

    TH1F* allHists[] = {
        sphericityHist_AllJets,
        sphericityHist_0Jets,
        sphericityHist_1Jet,
        sphericityHist_2Jets,
        sphericityHist_3PlusJets
    };

    for (auto hist : allHists) {
        for (int bin = 1; bin <= hist->GetNbinsX(); bin++) {
            double content = hist->GetBinContent(bin);
            if (content > 0) {  // Only consider positive values for log scale
                globalMax = max(globalMax, content);
                globalMin = min(globalMin, content);
            }
        }
    }

    // Add some padding to the ranges
    double yMin = (globalMin > 0) ? globalMin * 0.9 : 1e-1;  // Ensure positive min for log scale
    double yMax = globalMax * 1.1;  // Padding for upper limit

    // Set y-axis range for all histograms
    for (auto hist : allHists) {
        hist->GetYaxis()->SetRangeUser(yMin, yMax);
    }

    // Create individual canvases with consistent y-axis ranges
    for (auto hist : allHists) {
        TString canvasName = Form("c_%s", hist->GetName());
        TCanvas *canvas = new TCanvas(canvasName, canvasName, 800, 600);
        canvas->SetLogy();
        gPad->SetLeftMargin(0.15);
        hist->Draw("P E");
        hist->GetYaxis()->SetRangeUser(yMin, yMax);
    }

    // Save all plots to ROOT file
    TFile *outputFile = TFile::Open("output_sphericity.root", "RECREATE");
    if (!outputFile || outputFile->IsZombie()) {
        cout << "Error creating output file!" << endl;
        return;
    }
    
    // Write histograms to output file
    sphericityHist_AllJets->Write();
    sphericityHist_0Jets->Write();
    sphericityHist_1Jet->Write();
    sphericityHist_2Jets->Write();
    sphericityHist_3PlusJets->Write();
    rawJetMultiplicity->Write();
    filteredJetMultiplicity->Write();
    eventEtDist->Write();
    jetEtVsST->Write();
    
    // Close output file
    //outputFile->Close();
    
    cout << "Analysis complete!" << endl;
}

*/
