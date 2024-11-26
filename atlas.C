//this is the final code to produce the Aplanerity and the Sphericity from ATLAS paper

#include <iostream>
#include <algorithm>
#include <cmath>
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/PseudoJet.hh"

using namespace fastjet;
using namespace std;

// Define Matrix3x3 struct
struct Matrix3x3 {
    double m[3][3];
    
    Matrix3x3() {
        // Initialize the matrix to zero
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                m[i][j] = 0.0;
            }
        }
    }
    
    // Overload the + operator for matrix addition
    Matrix3x3 operator+(const Matrix3x3& other) const {
        Matrix3x3 result;
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                result.m[i][j] = m[i][j] + other.m[i][j];
            }
        }
        return result;
    }
};

// Function to calculate the determinant of a 2x2 matrix
double determinant2x2(double a, double b, double c, double d) {
    return a * d - b * c;
}

// Function to calculate the determinant of a 3x3 matrix
double determinant3x3(double matrix[3][3]) {
    double det = 0;
    for (int i = 0; i < 3; i++) {
        det += matrix[0][i] * determinant2x2(matrix[1][(i+1)%3], matrix[1][(i+2)%3], matrix[2][(i+1)%3], matrix[2][(i+2)%3]) * (i % 2 == 0 ? 1 : -1);
    }
    return det;
}

// Function to calculate the eigenvalues of a 3x3 matrix
void calculateEigenvalues(double matrix[3][3], double eigenvalues[3]) {
    double a = -1;
    double b = matrix[0][0] + matrix[1][1] + matrix[2][2];
    double c = determinant2x2(matrix[0][0], matrix[0][1], matrix[1][0], matrix[1][1]) +
               determinant2x2(matrix[0][0], matrix[0][2], matrix[2][0], matrix[2][2]) +
               determinant2x2(matrix[1][1], matrix[1][2], matrix[2][1], matrix[2][2]);
    double d = determinant3x3(matrix);

    // Solving the cubic equation a*x^3 + b*x^2 + c*x + d = 0
    double q = (3*c - b*b) / 9;
    double r = (9*b*c - 27*d - 2*b*b*b) / 54;
    double discriminant = q*q*q + r*r;

    if (discriminant > 0) {
        double s = cbrt(r + sqrt(discriminant));
        double t = cbrt(r - sqrt(discriminant));
        eigenvalues[0] = s + t - b / 3;
        eigenvalues[1] = eigenvalues[2] = 0; // Placeholder, as the full cubic solution is not implemented here.
    } else {
        double theta = acos(r / sqrt(-q*q*q));
        eigenvalues[0] = 2 * sqrt(-q) * cos(theta / 3) - b / 3;
        eigenvalues[1] = 2 * sqrt(-q) * cos((theta + 2*M_PI) / 3) - b / 3;
        eigenvalues[2] = 2 * sqrt(-q) * cos((theta + 4*M_PI) / 3) - b / 3;
    }
}

// Function to calculate the M_xyz matrix for a single jet
Matrix3x3 calculateMXYZ(const PseudoJet &jet) {
    Matrix3x3 M;
    double pX = jet.px();
    double pY = jet.py();
    double pZ = jet.pz();
    
    M.m[0][0] = pX * pX;
    M.m[0][1] = pX * pY;
    M.m[0][2] = pX * pZ;
    M.m[1][0] = pY * pX;
    M.m[1][1] = pY * pY;
    M.m[1][2] = pY * pZ;
    M.m[2][0] = pZ * pX;
    M.m[2][1] = pZ * pY;
    M.m[2][2] = pZ * pZ;
    
    return M;
}

void atlas() {
    // Open the ROOT file and get the tree
    TFile *inputFile = TFile::Open("/home/shalini/Software/pythia/examples/7_shalini.root");
    TTree *tree = (TTree *)inputFile->Get("Tree");
  
    // Set branch addresses to read particle data
    Int_t numParticles;
    Double_t px[500], py[500], pz[500], energy[500];
    tree->SetBranchAddress("numParticles", &numParticles);
    tree->SetBranchAddress("px", px);
    tree->SetBranchAddress("py", py);
    tree->SetBranchAddress("pz", pz);
    tree->SetBranchAddress("energy", energy);
  
    // Histograms for sphericity and aplanarity
    TH1F *sphericityHist = new TH1F("sphericityHist", "Sphericity;S;dN/dS", 100, 0, 5);
    TH1F *aplanarityHist = new TH1F("aplanarityHist", "Aplanarity;A;dN/dA", 100, 0, 5);
  
    // Jet clustering parameters
    double R = 0.4;
    JetDefinition jet_def(antikt_algorithm, R);
    vector<PseudoJet> particles;
  
    Int_t nentries = (Int_t)tree->GetEntries();
  
    for (Int_t i = 0; i < nentries; i++) {
        tree->GetEntry(i);
        particles.clear();
      
        // Prepare particles for jet clustering
        for (Int_t j = 0; j < numParticles; j++) {
            particles.push_back(PseudoJet(px[j], py[j], pz[j], energy[j]));
        }
      
        // Cluster jets
        ClusterSequence clustSeq(particles, jet_def);
        vector<PseudoJet> inclusive_jets = sorted_by_pt(clustSeq.inclusive_jets());
      
        // Ensure there are at least 3 jets
        if (inclusive_jets.size() < 3) continue;
      
        // Get the 1st, 2nd, and 3rd highest pt jets
        PseudoJet jet1 = inclusive_jets[0];
        PseudoJet jet2 = inclusive_jets[1];
        PseudoJet jet3 = inclusive_jets[2];
      
        // Calculate M_xyz matrices for these jets
        Matrix3x3 M1 = calculateMXYZ(jet1);
        Matrix3x3 M2 = calculateMXYZ(jet2);
        Matrix3x3 M3 = calculateMXYZ(jet3);
      
        // Sum the matrices
        Matrix3x3 Msum = M1 + M2 + M3;
      
        // Calculate the eigenvalues of the summed matrix
        double eigenvalues[3];
        calculateEigenvalues(Msum.m, eigenvalues);
      
        // Sort the eigenvalues in descending order
        sort(eigenvalues, eigenvalues + 3, greater<double>());
      
        // Normalize eigenvalues so that their sum is 1
        double sum = eigenvalues[0] + eigenvalues[1] + eigenvalues[2];
        eigenvalues[0] /= sum;
        eigenvalues[1] /= sum;
        eigenvalues[2] /= sum;
      
        // Calculate sphericity (S) and aplanarity (A)
        double S = 1.5 * (eigenvalues[1] + eigenvalues[2]);
        double A = 1.5 * eigenvalues[2];
      
        // Fill the histograms
        sphericityHist->Fill(S);
        aplanarityHist->Fill(A);
    }
  
    // Normalize histograms if needed
    sphericityHist->Scale(1.0 / sphericityHist->Integral());
    aplanarityHist->Scale(1.0 / aplanarityHist->Integral());

    //from here code is modified for the visibility of the plot
    
    // Sphericity Histogram Styles
    sphericityHist->SetMarkerStyle(20);      // Full circle
    sphericityHist->SetMarkerSize(1.2);      // Larger marker size
    sphericityHist->SetLineWidth(2);         // Thicker line
    sphericityHist->SetMarkerColor(kRed);    // Red markers
    sphericityHist->SetLineColor(kRed);      // Red lines

    // Aplanarity Histogram Styles
    aplanarityHist->SetMarkerStyle(21);      // Full square
    aplanarityHist->SetMarkerSize(1.2);      // Larger marker size
    aplanarityHist->SetLineWidth(2);         // Thicker line
    aplanarityHist->SetMarkerColor(kBlue);   // Blue markers
    aplanarityHist->SetLineColor(kBlue);     // Blue lines

    // **Create a Canvas for Plotting**
    TCanvas *c1 = new TCanvas("c1", "Sphericity and Aplanarity", 800, 600);
    c1->SetGrid(); // Add grid lines for better readability
    
    // **Draw the Histograms with Enhanced Styles**
    sphericityHist->Draw("E1"); // "E1" option draws error bars with markers
    aplanarityHist->Draw("E1 SAME"); // Overlay aplanarity histogram on the same canvas

    // **Add a Legend to Identify Histograms**
    TLegend *legend = new TLegend(0.7, 0.7, 0.9, 0.9); // Position (x1, y1, x2, y2)
    legend->SetBorderSize(0); // No border
    legend->SetFillStyle(0);   // Transparent fill
    legend->AddEntry(sphericityHist, "Sphericity", "lep"); // "lep" stands for Line, Error, Point
    legend->AddEntry(aplanarityHist, "Aplanarity", "lep");
    legend->Draw();
    
    c1->SaveAs("8_Sphericity_Aplanarity.png");
    
    //till here
    
    
    // Save histograms to a file
    TFile *outputFile = new TFile("atlaspaper_7mpion.root", "RECREATE");
    sphericityHist->Write();
    aplanarityHist->Write();
    //  outputFile->Close();
    
      return 0;
}
