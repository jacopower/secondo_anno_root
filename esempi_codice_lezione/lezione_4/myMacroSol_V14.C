void setStyle()
{
  gROOT->SetStyle("Plain");
  gStyle->SetPalette(57);
  gStyle->SetOptTitle(0);
}

Double_t myFunction(Double_t *x, Double_t *par)
{
  Double_t xx = x[0];
  Double_t val = par[0] * TMath::Exp(-(xx - par[1]) * (xx - par[1]) / 2. / par[2] / par[2]);
  return val;
}

void myMacro(Int_t nGen = 1E5)
{

  gStyle->SetOptStat(2210);
  gStyle->SetOptFit(1111);

  char *histName = new char[10];
  TH1F *h[2];
  for (Int_t i = 0; i < 2; i++)
  {
    sprintf(histName, "h%d", i);
    h[i] = new TH1F(histName, "test histogram", 100, -5, 5.);

    //cosmetics
    h[i]->SetMarkerStyle(20);
    h[i]->SetMarkerSize(0.5);
    h[i]->SetLineColor(1);
    h[i]->GetYaxis()->SetTitleOffset(1.2);
    h[i]->GetXaxis()->SetTitleSize(0.04);
    h[i]->GetYaxis()->SetTitleSize(0.04);
    h[i]->GetXaxis()->SetTitle("x");
    h[i]->GetYaxis()->SetTitle("Entries");
  }
  // filling histogram with predefined gaussian function

  h[0]->SetFillColor(4);
  h[1]->SetFillColor(2);

  h[0]->FillRandom("gaus", nGen); //gaus predefined function (G(0,1))

  TF1 *f = new TF1("f", myFunction, -10, 10, 3); //user defined function
  f->SetParameter(0, 1);
  f->SetParameter(1, 0);
  f->SetParameter(2, 1);
  h[1]->FillRandom("f", nGen);

  //Drawing

  TCanvas *c1 = new TCanvas("c1", "fitFunc examples", 200, 10, 600, 400);
  c1->Divide(1, 2);
  for (Int_t i = 0; i < 2; i++)
  {
    c1->cd(i + 1);
    h[i]->Draw("H");
    h[i]->Draw("E,P,SAME");
    c1->Print("testHisto.gif");
    c1->Print("testHisto.C");
    c1->Print("testHisto.root");
  }

  TCanvas *c2 = new TCanvas("c2", "fitFunc examples", 200, 10, 600, 400);
  h[0]->Fit("gaus");
  h[0]->Draw("H");
  h[0]->Draw("E,P,SAME");

  TCanvas *c3 = new TCanvas("c3", "fitFunc examples", 200, 10, 600, 400);
  f->SetParameter(0, 4000);
  f->SetParameter(1, 0);
  f->SetParameter(2, 1);
  f->SetLineColor(kCyan);
  TFitResultPtr fRes = h[1]->Fit(f, "S");
  h[1]->Draw("H");
  h[1]->Draw("E,P,SAME");
  TLegend *leg = new TLegend(.1, .7, .3, .9, "test Fit ");
  leg->SetFillColor(0);
  leg->AddEntry(h[1], "Punti sperimentali");
  leg->AddEntry(f, "Fit Gaussiano");
  leg->Draw("S");

  h[1]->GetListOfFunctions()->Print();
  TF1 *fitFunc = h[1]->GetFunction("f");

  Double_t mean = fitFunc->GetParameter(1);
  Double_t meanErr = fitFunc->GetParError(1);
  Double_t sigma = fitFunc->GetParameter(2);
  Double_t sigmaErr = fitFunc->GetParError(2);
  Double_t chiSquare = fitFunc->GetChisquare();
  Double_t nDOF = fitFunc->GetNDF();
  Double_t Prob = fitFunc->GetProb();
  cout << "Mean = " << mean << " +/- " << meanErr << endl;
  cout << "Sigma = " << sigma << " +/- " << sigmaErr << endl;
  cout << "ChiSquare = " << chiSquare << " , nDOF " << nDOF << endl;
  cout << "ChiSquare Probability= " << Prob << endl;

  //Covariance and Correlation

  TMatrixD cov = fRes->GetCovarianceMatrix();
  TMatrixD cor = fRes->GetCorrelationMatrix();
  cov.Print();
  cor.Print();
}
