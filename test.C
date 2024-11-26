/*
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/PseudoJet.hh"
#include "fastjet/JetDefinition.hh"
#include "fastjet/Selector.hh"

using namespace std;
using namespace fastjet;

double calculateTransverseSphericity(const vector<double> &px, const vector<double> &py) {
    double sxx = 0, syy = 0, sxy = 0;
    double sum_pt = 0;
    for (size_t i = 0; i < px.size(); ++i) {
        double pt = sqrt(px[i] * px[i] + py[i] * py[i]);
        sum_pt += pt;
        if (pt > 0) {
            sxx += px[i] * px[i] / pt;
            syy += py[i] * py[i] / pt;
            sxy += px[i] * py[i] / pt;
        }
    }

    if (sum_pt == 0) return -1;  // Invalid case: zero sum_pt leads to division by zero

    // Normalize by sum_pt
    sxx /= sum_pt;
    syy /= sum_pt;
    sxy /= sum_pt;

    double trace = sxx + syy;
    double det = sxx * syy - sxy * sxy;

    double lambda1 = (trace + sqrt(trace * trace - 4 * det)) / 2;
    double lambda2 = (trace - sqrt(trace * trace - 4 * det)) / 2;

    return 2 * lambda2 / (lambda1 + lambda2);  // Transverse Sphericity S_T
}

void test() {
    // Open the ROOT file and get the tree
    TFile *inputFile = TFile::Open("/home/shalini/Software/pythia/examples/13.6mpion_shalini.root");
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
    int zeroPtEvents = 0;

    for (Int_t i = 0; i < nentries; i++) {
        tree->GetEntry(i);
        particles.clear();
        vector<double> event_px, event_py;

        bool hasNonZeroPt = false;

        // Prepare particles for jet clustering and sphericity calculation
        for (Int_t j = 0; j < numParticles; j++) {
            double pt_value = sqrt(px[j] * px[j] + py[j] * py[j]);
            if (pt_value > 0) hasNonZeroPt = true;  // Check for non-zero pt
            
            particles.push_back(PseudoJet(px[j], py[j], pz[j], energy[j]));
            event_px.push_back(px[j]);
            event_py.push_back(py[j]);
        }

        // If no particle has non-zero pt, skip the event
        if (!hasNonZeroPt) {
            zeroPtEvents++;
            logFile << "Skipping event " << i << " with all zero pt particles." << endl;
            continue;
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
            sphericityHist_AllJets->Fill(ST);

            if (numJets == 0) sphericityHist_0Jets->Fill(ST);
            else if (numJets == 1) sphericityHist_1Jet->Fill(ST);
            else if (numJets == 2) sphericityHist_2Jets->Fill(ST);
            else if (numJets >= 3) sphericityHist_3PlusJets->Fill(ST);
        }
    }

    // Print the number of skipped zero pt events
    cout << "Number of events skipped due to zero pt: " << zeroPtEvents << endl;

    // Close the log file
    logFile.close();

    // Normalize histograms
    double scale = 1.0 / validEntries;
    sphericityHist_AllJets->Scale(scale);
    sphericityHist_0Jets->Scale(scale);
    sphericityHist_1Jet->Scale(scale);
    sphericityHist_2Jets->Scale(scale);
    sphericityHist_3PlusJets->Scale(scale);

    // Clean up and save histograms
    TFile *outputFile = new TFile("output_sphericity.root", "RECREATE");
    sphericityHist_AllJets->Write();
    sphericityHist_0Jets->Write();
    sphericityHist_1Jet->Write();
    sphericityHist_2Jets->Write();
    sphericityHist_3PlusJets->Write();
    //outputFile->Close();

    //delete inputFile;
}


*/

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/PseudoJet.hh"
#include "fastjet/JetDefinition.hh"
#include "fastjet/Selector.hh"

using namespace std;
using namespace fastjet;

double calculateTransverseSphericity(const vector<double> &px, const vector<double> &py) {
    double sxx = 0, syy = 0, sxy = 0;
    double sum_pt = 0;
    for (size_t i = 0; i < px.size(); ++i) {
        double pt = sqrt(px[i] * px[i] + py[i] * py[i]);
        sum_pt += pt;
        if (pt > 0) {
            sxx += px[i] * px[i] / pt;
            syy += py[i] * py[i] / pt;
            sxy += px[i] * py[i] / pt;
        }
    }

    if (sum_pt == 0) return -1;  // Invalid case: zero sum_pt leads to division by zero

    // Normalize by sum_pt
    sxx /= sum_pt;
    syy /= sum_pt;
    sxy /= sum_pt;

    double trace = sxx + syy;
    double det = sxx * syy - sxy * sxy;

    double lambda1 = (trace + sqrt(trace * trace - 4 * det)) / 2;
    double lambda2 = (trace - sqrt(trace * trace - 4 * det)) / 2;

    return 2 * lambda2 / (lambda1 + lambda2);  // Transverse Sphericity S_T
}

void test() {
    // Open the ROOT file and get the tree
    TFile *inputFile = TFile::Open("/home/shalini/Software/pythia/examples/13.6mpion_shalini.root");
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

    // Histogram for leading hadron pt
    TH1F *leadingPtHist = new TH1F("leadingPtHist", "Leading Hadron p_{T};p_{T} (GeV);dN/dp_{T}", 100, 0, 100); // Adjust the binning as needed

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
    int zeroPtEvents = 0;

    for (Int_t i = 0; i < nentries; i++) {
        tree->GetEntry(i);
        particles.clear();
        vector<double> event_px, event_py;

        bool hasNonZeroPt = false;

        // Prepare particles for jet clustering and sphericity calculation
        for (Int_t j = 0; j < numParticles; j++) {
            double pt_value = sqrt(px[j] * px[j] + py[j] * py[j]);
            if (pt_value > 0) hasNonZeroPt = true;  // Check for non-zero pt
            
            particles.push_back(PseudoJet(px[j], py[j], pz[j], energy[j]));
            event_px.push_back(px[j]);
            event_py.push_back(py[j]);
        }

        // If no particle has non-zero pt, skip the event
        if (!hasNonZeroPt) {
            zeroPtEvents++;
            logFile << "Skipping event " << i << " with all zero pt particles." << endl;
            continue;
        }

        // Find leading hadron pt
        double leadingPt = 0;
        for (Int_t j = 0; j < numParticles; j++) {
            double pt_value = sqrt(px[j] * px[j] + py[j] * py[j]);
            if (pt_value > leadingPt) {
                leadingPt = pt_value;
            }
        }
        leadingPtHist->Fill(leadingPt); // Fill the leading hadron pt histogram

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
            sphericityHist_AllJets->Fill(ST);

            if (numJets < 1) sphericityHist_0Jets->Fill(ST);
            else if (numJets == 1) sphericityHist_1Jet->Fill(ST);
            else if (numJets == 2) sphericityHist_2Jets->Fill(ST);
            else if (numJets >= 3) sphericityHist_3PlusJets->Fill(ST);
        }
    }

    // Print the number of skipped zero pt events
    cout << "Number of events skipped due to zero pt: " << zeroPtEvents << endl;

    // Close the log file
    logFile.close();

    // Normalize histograms
    double scale = 1.0 / validEntries;
    sphericityHist_AllJets->Scale(scale);
    sphericityHist_0Jets->Scale(scale);
    sphericityHist_1Jet->Scale(scale);
    sphericityHist_2Jets->Scale(scale);
    sphericityHist_3PlusJets->Scale(scale);

    // Clean up and save histograms
    TFile *outputFile = new TFile("output_sphericity.root", "RECREATE");
    sphericityHist_AllJets->Write();
    sphericityHist_0Jets->Write();
    sphericityHist_1Jet->Write();
    sphericityHist_2Jets->Write();
    sphericityHist_3PlusJets->Write();
    leadingPtHist->Write(); // Save the leading hadron pt histogram
    //outputFile->Close();

    //delete inputFile;
}
