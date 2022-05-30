void prova()
{
  TCanvas *c1 = new TCanvas("c1", "The FillRandom example",
                            200, 10, 700, 900);
  c1->SetFillColor(18);
  TPad *pad1 = new TPad("pad1", "The pad with the function",
                       0.05, 0.50, 0.95, 0.95, 21);
  TPad *pad2 = new TPad("pad2", "The pad with the histogram",
                       0.05, 0.05, 0.95, 0.45, 21);
  pad1->Draw();
  pad2->Draw();
}