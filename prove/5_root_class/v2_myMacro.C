{
  TF1 *fIn = new TF1("fIn", "x", 0., 1.);
  TH1F *hIn = new TH1F("hIn", "input histo", 100, 0., 1.);

  myClass A(hIn);

  A.Generate(fIn, 10000);
  A.ShowHisto();

  TFile *file = new TFile("test.root", "RECREATE");
  A.Write("A");
  file->Close();
}