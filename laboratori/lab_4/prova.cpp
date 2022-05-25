#include <iostream>
#include <vector>
#include "TF1.h"
#include "TMath.h"
#include "TGraph.h"
#include "TCanvas.h"

using namespace std;
void Solve(vector<double> &roots, int nest, TF1 *f, double s, double e)
{
  double r = f->GetX(0, s, e);
  cout << "(" << nest << ") start=" << s << " end=" << e << " root=" << r << endl;
  roots.push_back(r);
  if (r != s && r != e)
  {
    Solve(roots, nest + 1, f, s, r);
    Solve(roots, nest + 1, f, r, e);
  }
}
void FindRoot()
{
  double s = -TMath::PiOver2();
  double e = 2 * TMath::TwoPi() + TMath::PiOver2();
  cout << "start =" << s << endl;
  cout << "end   =" << e << endl;
  // Create the function and wrap it
  TF1 *f = new TF1("Sin Function", "sin(x)", s, e);
  f->Draw();
  vector<double> roots;
  Solve(roots, 0, f, s, e);
  sort(roots.begin(), roots.end());
  double *xx = new double[roots.size()];
  double *yy = new double[roots.size()];
  // If this is put here; then only the first point
  // is shown in the plot.
  // TGraph* gr = new TGraph(7, xx, yy);
  cout << "Found " << roots.size() << " roots:";
  for (int i = 0; i < roots.size(); ++i)
  {
    xx[i] = roots[i];
    yy[i] = 0;
    cout << xx[i] << " ";
  }
  cout << endl;
  TGraph *gr = new TGraph(7, xx, yy);
  gr->SetMarkerColor(kRed);
  gr->SetMarkerStyle(21);
  gr->Draw("P");
  c1->Update();
  cout.flush();
}