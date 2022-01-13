void setStyle()
{
  gROOT->SetStyle("Plain");
  gStyle->SetPalette(57);
  gStyle->SetOptTitle(0);
}

void myMacro()
{

  TH1F *h = new TH1F("h", "test histogram", 100, -5., 5.);

  // filling histogram with predefined gaussian function
  h->FillRandom("gaus", 1E6);

  h->Draw();
}
