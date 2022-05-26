
TCanvas *pavetext()
{
  TCanvas *c = new TCanvas("c");
  TPaveText *pt = new TPaveText(.05, .1, .95, .8);

  pt->AddText("A TPaveText can contain severals line of text.");
  pt->AddText("They are added to the pave using the AddText method.");
  pt->AddLine(.0, .5, 1., .5);
  pt->AddText("Even complex TLatex formulas can be added:");
  TText *t1 = pt->AddText("F(t) = #sum_{i=-#infty}^{#infty}A(i)cos#[]{#frac{i}{t+i}}");

  t1->SetTextColor(kBlue);

  pt->Draw();

  TText *t2 = pt->GetLineWith("Even");
  t2->SetTextColor(kOrange + 1);

  return c;
}

//**************************************************************
TCanvas *statsEditing() {
   // Create and plot a test histogram with stats
   TCanvas *se = new TCanvas;
   TH1F *h = new TH1F("h","test",100,-3,3);
   h->FillRandom("gaus",3000);
   gStyle->SetOptStat();
   h->Draw();
   se->Update();
 
   // Retrieve the stat box
   TPaveStats *ps = (TPaveStats*)se->GetPrimitive("stats");
   ps->SetName("mystats");
   TList *listOfLines = ps->GetListOfLines();
 
   // Remove the RMS line
   TText *tconst = ps->GetLineWith("RMS");
   listOfLines->Remove(tconst);
 
   // Add a new line in the stat box.
   // Note that "=" is a control character
   TLatex *myt = new TLatex(0,0,"Test = 10");
   myt ->SetTextFont(42);
   myt ->SetTextSize(0.04);
   myt ->SetTextColor(kRed);
   listOfLines->Add(myt);
 
   // the following line is needed to avoid that the automatic redrawing of stats
   h->SetStats(0);
 
   se->Modified();
   return se;
}
