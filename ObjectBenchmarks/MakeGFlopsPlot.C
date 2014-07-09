#include <iostream>
#include <fstream>
#include <vector>
#include <utility>
#include <string>
#include "TGraph.h"
#include "TCanvas.h"
#include "TLegend.h"


using namespace std;

//vector<pair<int,float> > getTiming(const string& infilename);

//TGraph *makeGraph(const vector<pair<int,float> > timingvec, const string& name);
TGraph *makeGraph(const string& filename, const string& title, const EColor& color, const EMarkerStyle& shape);
TCanvas *drawGraphs(vector<TGraph*> graphs);

int main(){
  
  //  vector<pair<int,float> > test = getTiming("output/out_2arrays_4calc.txt");



  vector<TGraph*> graphs_host;
  graphs_host.push_back(makeGraph("output/out_2arrays_5calc.txt", "2 Arrays, 5 Calculations", kRed, kCircle));
  graphs_host.push_back(makeGraph("output/out_3arrays_5calc.txt", "3 Arrays, 5 Calculations", kBlue, kOpenTriangleUp));
  graphs_host.push_back(makeGraph("output/out_4arrays_5calc.txt", "4 Arrays, 5 Calculations", kCyan, kFullSquare));
  TCanvas *narrays_5calcs_host=drawGraphs(graphs_host);


  vector<TGraph*> graphs_mic;
  graphs_mic.push_back(makeGraph("output_mic/out_2arrays_5calc.txt", "2 Arrays, 5 Calculations", kRed, kCircle));
  graphs_mic.push_back(makeGraph("output_mic/out_3arrays_5calc.txt", "3 Arrays, 5 Calculations", kBlue, kOpenTriangleUp));
  graphs_mic.push_back(makeGraph("output_mic/out_4arrays_5calc.txt", "4 Arrays, 5 Calculations", kCyan, kFullSquare));
  TCanvas *narrays_5calcs_mic=drawGraphs(graphs_mic);
  
  

  vector<TGraph*> resid_host;
  resid_host.push_back(makeGraph("output/resid_6arrays_12calc.txt", "6 Arrays, 12 Calculations", kRed, kCircle));
  resid_host.push_back(makeGraph("output/resid_s6floats_12calc.txt", "Array of SoF, 12 Calculations", kBlue, kOpenTriangleUp));
  resid_host.push_back(makeGraph("output/resid_shsa_12calc.txt", "Array of SoSA, 12 Calculations", kCyan, kFullSquare));
  TCanvas *resid_12calcs_host=drawGraphs(resid_host);

  vector<TGraph*> resid_mic;
  resid_mic.push_back(makeGraph("output_mic/resid_6arrays_12calc.txt", "6 Arrays, 12 Calculations", kRed, kCircle));
  resid_mic.push_back(makeGraph("output_mic/resid_s6floats_12calc.txt", "Array of SoF, 12 Calculations", kBlue, kOpenTriangleUp));
  resid_mic.push_back(makeGraph("output_mic/resid_shsa_12calc.txt", "Array of SoSA, 12 Calculations", kCyan, kFullSquare));
  TCanvas *resid_12calcs_mic=drawGraphs(resid_mic);

  return 0;
}




// vector<pair<int,float> > getTiming(const string& infilename){
//   ifstream inf;
//   inf.open(infilename.c_str());

//   vector<pair<int,float> > inv;

//   while(inf.good()) {
// 	int tempi=0;
// 	float tempf=0;
// 	inf>>tempi>>tempf;
// 	if(inf.eof()){ break; }
// 	cout<<tempi<<" "<<tempf<<endl;
// 	inv.push_back(pair<int,float>(tempi, tempf));
//   }

//   return inv;
// }

TGraph *makeGraph(const string& filename, const string& title, const EColor& color, const EMarkerStyle& shape){
  TGraph *temp = new TGraph(filename.c_str());
  //  temp->SetDirectory(NULL);
  temp->SetTitle(title.c_str());
  //temp->SetMarkerColor(color);
  temp->SetLineColor(color);
  //  temp->SetMarkerStyle(shape);
  return temp;
}

//TGraph *makeGraph(const vector<pair<int,float> > timingvec, const string& name){


double TGraphMaxY(TGraph *g){
  double maximum = 0;
  double *y = g->GetY();
  for(int i=0; i<g->GetN(); ++i){
	maximum = std::max(maximum, y[i]);
	//	cout<<y[i]<<endl;
  }
  //  cout<<"graph max "<<maximum<<endl;
  return maximum;
}

double TGraphsMaxY(vector<TGraph*> graphs){
  double maximum = 0;
  for(unsigned int i=0; i<graphs.size(); ++i){
	//	cout<<"graph "<<i<<endl;
	maximum = std::max(maximum, TGraphMaxY(graphs[i]));
  }
  //  cout<<"graphs max "<<maximum<<endl;
  return maximum;
}



TCanvas *drawGraphs(vector<TGraph*> graphs){
  TCanvas *c = new TCanvas();
  c->cd();
  c->SetLogx();

  TLegend *l = new TLegend(0.65,0.7,0.9,0.9);
  double maximum = TGraphsMaxY(graphs);

  if(graphs.size()>=1){ 
	graphs[0]->Draw();
	graphs[0]->GetYaxis()->SetRangeUser(0., maximum*1.1);
	graphs[0]->Draw();
	l->AddEntry(graphs[0], graphs[0]->GetTitle(), "l"); 
  }
  for(unsigned int i=1; i<graphs.size(); ++i){
	graphs[i]->Draw("same");
	l->AddEntry(graphs[i],graphs[i]->GetTitle(), "l");
  }
  l->Draw();
  return c;
}



//}
