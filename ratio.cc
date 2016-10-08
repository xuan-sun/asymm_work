#include 	 "BetaSpectrum.hh"

#include	 <iostream>
#include	 <fstream>
#include	 <TGaxis.h>
#include	 <sstream>
#include	 <TGraph.h>
#include	 <TGraphErrors.h>
#include	 <TCanvas.h>
#include	 <TApplication.h>
#include	 <stdlib.h>
#include	 <TF1.h>
#include	 <TH1.h>
#include	 <TProfile.h>
#include	 <TObjArray.h>
#include	 <TStyle.h>
#include	 <TMarker.h>
#include	 <math.h>
#include	 <TStyle.h>
#include	 <TPaveStats.h>
#include	 <TPaveText.h>
#include	 <vector>
#include	 <string.h>
#include	 <fstream>
#include	 <TROOT.h>
#include	 <TFile.h>
#include	 <TLegend.h>
#include         <TLegendEntry.h>
#include	 <time.h>
#include	 <TH2F.h>
#include         <assert.h>
#include	 <string>
#include	 <TRandom.h>
#include 	 <TTree.h>
#include	 <TChain.h>
#include	 <TVector.h>
using		 namespace std;

//required later for plot_program
TApplication plot_program("FADC_readin",0,0,0,0);

void FillTheoryCorrections(vector <double> *Eth, vector <double> *Ath, vector <double> *Eth_error, vector <double> *Ath_error);
void LoadFile(TString fileName, vector <double> *E, vector <double> *A, vector <double> *E_error, vector <double> *A_error);
void PlotHist(TCanvas *C, int styleIndex, int canvaxIndex, TH1D *hPlot, TString title, TString command);
void PlotGraph(TCanvas *C, int styleIndex, int canvasIndex, TGraphErrors* gPlot, TString title, TString command);

int main(int argc, char* argv[])
{
  if(argc < 4)
  {
    cout << "Error: improper input. Must give:" << endl;
    cout << "(executable) (Asymm Octet.txt Path) (octet #) (analysis choice)" << endl;
    return 0;
  }

  // read in arguments. Pass it a /Data/ folder and number of TChains to store
  Int_t octNb = atoi(argv[2]);
  Int_t analysis = atoi(argv[3]);
  TString path = argv[1];

  TCanvas *C = new TCanvas("canvas", "canvas");
//  C -> Divide(2,1);
  gROOT -> SetStyle("Plain");   //on my computer this sets background to white, finally!

//  TString fileName1 = TString::Format("%sAsymm_Octet-%i_Analysis-%i_13.txt", path.Data(), octNb, analysis);
  TString fileName2 = TString::Format("%sAsymm_Octet-%i_Analysis-%i_20.txt", path.Data(), octNb, analysis);

  vector <double> *A_pure = new vector <double>;
  vector <double> *A_error_pure = new vector <double>;
  vector <double> *E_pure = new vector <double>;
  vector <double> *E_error_pure = new vector<double>;

  vector <double> *A_proc = new vector <double>;
  vector <double> *A_error_proc = new vector <double>;
  vector <double> *E_proc = new vector <double>;
  vector <double> *E_error_proc = new vector<double>;

//  LoadFile(fileName1, E_pure, A_pure, E_error_pure, A_error_pure);
  FillTheoryCorrections(E_pure, A_pure, E_error_pure, A_error_pure);
  LoadFile(fileName2, E_proc, A_proc, E_error_proc, A_error_proc);


  vector <double> pureOverProc_ratio;
  vector <double> ratioError;
  double ratio = 0;
  double ratioE = 0;
  for(unsigned int i = 0; i < A_pure->size(); i++)
  {
    if((*A_pure)[i] == 0 || (*A_proc)[i] == 0)
    {
      ratio = 0;
    }
    else
    {
      ratio = (-1)*(*A_pure)[i]/ (*A_proc)[i];
    }

    if(ratio == 0 || (*A_pure)[i] == 0 || (*A_proc)[i] == 0)
    {
      ratioE = 0;
    }
    else
    {
      // error when there is statistical error in each pure and proc.
      // When error in pure is 0 though it will still propagate correctly.
      ratioE = abs(ratio)*sqrt( ((*A_error_pure)[i] / (*A_pure)[i])*((*A_error_pure)[i] / (*A_pure)[i]) +
		((*A_error_proc)[i] / (*A_proc)[i])*((*A_error_proc)[i] / (*A_proc)[i]));
    }

    pureOverProc_ratio.push_back(ratio);
    ratioError.push_back(ratioE);
  }
  //	Creates output .txt files containing the correction ratios
  ofstream outfile;
  char tmp[250];
  sprintf(tmp, "ThOverProc_Octet-%i_Analysis-%i.txt", octNb, analysis);
  outfile.open(tmp, ios::app);
  outfile << "Energy midpoint \t Ratio (pure/proc) \t E error \t Ratio error" << endl;
  for(unsigned h = 0; h < pureOverProc_ratio.size(); h++)
  {
    // any index from 0 to 7 will do, since they all are the same energy bins
    outfile << (*E_pure)[h] << "\t"
            << pureOverProc_ratio[h] << "\t"
            << (*E_error_pure)[h] << "\t"
            << ratioError[h] << "\n";
  }
  outfile.close();

/*  // makes a TGraphErrors and plots it/saves to pdf/allows instant viewing
  TVectorD xValue(E_pure->size(), &(*E_pure)[0]);
  TVectorD yValue(pureOverProc_ratio.size(), &(pureOverProc_ratio)[0]);
  TVectorD xErrors(E_error_pure->size(), &(*E_error_pure)[0]);
  TVectorD yErrors(ratioError.size(), &(ratioError)[0]);
  TGraphErrors* graphA = new TGraphErrors(xValue, yValue, xErrors, yErrors);
  PlotGraph(C, 1, 1, graphA, TString::Format("Asymm ratio for Octet-%i, Analysis-%i", octNb, analysis), "AP");
  C -> Print(TString::Format("RatioPlot_Octet-%i_Analysis-%i.pdf", octNb, analysis));
  plot_program.Run();
*/
/*  // Creates MPM-like ratio plots and either prints to a .pdf or saves to a .root Tfile.
  double sum = 0;
  double sumE = 0;
  vector <double> yVal, yErr, xVal, xErr;
  TString name = TString::Format("Asymm_ratio_Octet-%i_Analysis-%i.root", octNb, analysis);
  TFile f(name, "RECREATE");
  for(int t = 0; t < pureOverProc_ratio.size(); t++)
  {
    sum = sum + pureOverProc_ratio[t];
    sumE = sumE + ratioError[t];

    if(t%6 == 0)
    {
      yVal.push_back(sum/6.);
      yErr.push_back(sumE/6.);
      xVal.push_back(t*10 - 30);
      xErr.push_back(5);

      sum = 0;
      sumE = 0;
    }
  }
  TVectorD xValue(xVal.size(), &(xVal)[0]);
  TVectorD yValue(yVal.size(), &(yVal)[0]);
  TVectorD xErrors(xErr.size(), &(xErr)[0]);
  TVectorD yErrors(yErr.size(), &(yErr)[0]);
  TGraphErrors* graphA = new TGraphErrors(xValue, yValue, xErrors, yErrors);
//  PlotGraph(C, 1, 1, graphA, TString::Format("Asymm ratio for Octet-%i, Analysis-%i", octNb, analysis), "AXL");
//  C -> Print(TString::Format("MPM_bins_RatioPlot_Octet-%i_Analysis-%i.pdf", octNb, analysis));
//  plot_program.Run();
  graphA -> Write(name);
*/

  return 0;
}

void FillTheoryCorrections(vector <double> *Eth, vector <double> *Ath, vector <double> *Eth_error, vector <double> *Ath_error)
{
  for(int e = 5; e < 1201; e = e + 10)
  {
    Eth->push_back(e);
//    A->push_back(correctedAsymmetry(e, 0.5)*beta(e)/2.);
    Ath->push_back(A0_PDG*asymmetryCorrectionFactor(e)*beta(e)/2.);
    Eth_error->push_back(5);
    Ath_error->push_back(0);
  }

}


void LoadFile(TString fileName, vector <double> *E, vector <double> *A, vector <double> *E_error, vector <double> *A_error)
{
  string buf;
  ifstream infile;
  cout << "The file being opened is: " << fileName << endl;
  infile.open(fileName);

  if(!infile.is_open())
    cout << "Problem opening " << fileName << endl;

  double a = 0;
  double ae = 0;
  double e = 0;
  double ee = 0;
  int lineNb = 0;

  while(getline(infile, buf))
  {
    lineNb++;

    istringstream bufstream(buf);
    if(!bufstream.eof())
    {
      bufstream >> e >> a >> ee >> ae;

      if(lineNb == 1)
      {
        continue;
      }

      E->push_back(e);
      A->push_back(a);
      E_error->push_back(ee);
      A_error->push_back(ae);
    }
  }
}

void PlotHist(TCanvas *C, int styleIndex, int canvasIndex, TH1D *hPlot, TString title, TString command)
{
  C -> cd(canvasIndex);
  hPlot -> SetTitle(title);
  hPlot -> GetXaxis() -> SetTitle("Energy (keV)");
  hPlot -> GetXaxis() -> CenterTitle();
  hPlot -> GetYaxis() -> SetTitle("Hits Rate/Ratio (N/Total)");
  hPlot -> GetYaxis() -> CenterTitle();
//  hPlot -> GetYaxis() -> SetRangeUser(-0.15, 0.15);

  if(styleIndex == 1)
  {
    hPlot -> SetFillColor(46);
    hPlot -> SetFillStyle(3004);
  }
  if(styleIndex == 2)
  {
    hPlot -> SetFillColor(38);
    hPlot -> SetFillStyle(3005);
  }

  hPlot -> Draw(command);
}

void PlotGraph(TCanvas *C, int styleIndex, int canvasIndex, TGraphErrors* gPlot, TString title, TString command)
{
  C -> cd(canvasIndex);
  gPlot -> SetTitle(title);
  gPlot -> GetXaxis() -> SetTitle("Energy (keV)");
  gPlot -> GetXaxis() -> SetRangeUser(0, 800);
  gPlot -> GetXaxis() -> CenterTitle();
  gPlot -> GetYaxis() -> SetTitle("Correction to A");
  gPlot -> GetYaxis() -> CenterTitle();
  gPlot -> SetMinimum(-0.05);
  gPlot -> SetMaximum(0.05);

  if(styleIndex == 1)
  {
    gPlot -> SetMarkerColor(46);
  }
  if(styleIndex == 2)
  {
    gPlot -> SetMarkerColor(38);
  }
  if(styleIndex == 3)
  {
    gPlot -> SetMarkerColor(30);
  }
  gPlot -> SetMarkerStyle(7);
  gPlot -> SetMarkerSize(1);

//  TLine *line = new TLine(0, 1, 1200, 1);

  gPlot -> Draw(command);
//  line -> SetLineStyle(7);
//  line -> Draw("SAME");

  TLine *line2 = new TLine(0, 0, 800, 0);
  line2 -> Draw("SAME");

}


