
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