void c1()
{
//=========Macro generated from canvas: c1/c1
//=========  (Thu May 26 19:40:31 2022) by ROOT version 6.22/08
   TCanvas *c1 = new TCanvas("c1", "c1",2160,75,1745,1009);
   c1->ToggleEventStatus();
   c1->ToggleToolTips();
   c1->Range(0,0,1,1);
   c1->SetFillColor(0);
   c1->SetBorderMode(0);
   c1->SetBorderSize(2);
   c1->SetFrameBorderMode(0);
   
   TPaveText *pt = new TPaveText(0.1162669,0.6416244,0.4734098,0.9096447,"br");
   pt->SetLabel("ddd");
   pt->SetBorderSize(2);
   pt->SetFillColor(38);
   pt->SetTextColor(2);
   TLine *pt_Line = pt->AddLine(0,-2.393939,0,-2.393939);
   pt_Line = pt->AddLine(0,-2.393939,0,-2.393939);
   TText *pt_LaTex = pt->AddText("ciao allor asono dafa");
   pt_LaTex = pt->AddText("r =");
   pt->Draw();
   c1->Modified();
   c1->cd();
   c1->SetSelected(c1);
   c1->ToggleToolBar();
}
