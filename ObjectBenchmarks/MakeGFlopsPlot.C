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



  // vector<TGraph*> graphs_host;
  // graphs_host.push_back(makeGraph("output/out_2arrays_5calc.txt", "2 Arrays, 5 Calculations", kRed, kCircle));
  // graphs_host.push_back(makeGraph("output/out_3arrays_5calc.txt", "3 Arrays, 5 Calculations", kBlue, kOpenTriangleUp));
  // graphs_host.push_back(makeGraph("output/out_4arrays_5calc.txt", "4 Arrays, 5 Calculations", kCyan, kFullSquare));
  // TCanvas *narrays_5calcs_host=drawGraphs(graphs_host);


  // vector<TGraph*> graphs_mic;
  // graphs_mic.push_back(makeGraph("output_mic/out_2arrays_5calc.txt", "2 Arrays, 5 Calculations", kRed, kCircle));
  // graphs_mic.push_back(makeGraph("output_mic/out_3arrays_5calc.txt", "3 Arrays, 5 Calculations", kBlue, kOpenTriangleUp));
  // graphs_mic.push_back(makeGraph("output_mic/out_4arrays_5calc.txt", "4 Arrays, 5 Calculations", kCyan, kFullSquare));
  // TCanvas *narrays_5calcs_mic=drawGraphs(graphs_mic);
  
  

  // vector<TGraph*> resid_host;
  // resid_host.push_back(makeGraph("output/resid_6arrays_12calc.txt", "6 Arrays, 12 Calculations", kRed, kCircle));
  // resid_host.push_back(makeGraph("output/resid_s6floats_12calc.txt", "Array of SoF, 12 Calculations", kBlue, kOpenTriangleUp));
  // resid_host.push_back(makeGraph("output/resid_shsa_12calc.txt", "Array of SoSA, 12 Calculations", kCyan, kFullSquare));
  // TCanvas *resid_12calcs_host=drawGraphs(resid_host);

  // vector<TGraph*> resid_mic;
  // resid_mic.push_back(makeGraph("output_mic/resid_6arrays_12calc.txt", "6 Arrays, 12 Calculations", kRed, kCircle));
  // resid_mic.push_back(makeGraph("output_mic/resid_s6floats_12calc.txt", "Array of SoF, 12 Calculations", kBlue, kOpenTriangleUp));
  // resid_mic.push_back(makeGraph("output_mic/resid_shsa_12calc.txt", "Array of SoSA, 12 Calculations", kCyan, kFullSquare));
  // TCanvas *resid_12calcs_mic=drawGraphs(resid_mic);


  vector<TGraph*> resid_dl_2097152;
  vector<TGraph*> resid_dl_262144;
  vector<TGraph*> resid_dl_32768;
  vector<TGraph*> resid_dl_4096;
  vector<TGraph*> resid_dl_512;
  vector<TGraph*> resid_dl_64;
  //  vector<TGraph*> resid_dl_8;

  resid_dl_2097152.push_back(makeGraph("output/resid_dl_shsa_e2097152.txt", "host shsa 2097152", kRed, kCircle));
  resid_dl_262144.push_back(makeGraph("output/resid_dl_shsa_e262144.txt", "host shsa 262144", kRed, kCircle));
  resid_dl_32768.push_back(makeGraph("output/resid_dl_shsa_e32768.txt", "host shsa 32768", kRed, kCircle));
  resid_dl_4096.push_back(makeGraph("output/resid_dl_shsa_e4096.txt", "host shsa 4096", kRed, kCircle)); 
  resid_dl_512.push_back(makeGraph("output/resid_dl_shsa_e512.txt", "host shsa 512", kRed, kCircle));
  resid_dl_64.push_back(makeGraph("output/resid_dl_shsa_e64.txt", "host shsa 64", kRed, kCircle));
  //  resid_dl_8.push_back(makeGraph("output/resid_dl_shsa_e8.txt", "host shsa 8", kRed, kCircle));
  
  resid_dl_2097152.push_back(makeGraph("output/resid_dl_6arrays_e2097152.txt", "host a 2098152", kBlue, kOpenTriangleDown));
  resid_dl_262144.push_back(makeGraph("output/resid_dl_6arrays_e262144.txt", "host a 262144", kBlue, kOpenTriangleDown));
  resid_dl_32768.push_back(makeGraph("output/resid_dl_6arrays_e32768.txt", "host a 32768", kBlue, kOpenTriangleDown));
  resid_dl_4096.push_back(makeGraph("output/resid_dl_6arrays_e4096.txt", "host a 4096", kBlue, kOpenTriangleDown)); 
  resid_dl_512.push_back(makeGraph("output/resid_dl_6arrays_e512.txt", "host a 512", kBlue, kOpenTriangleDown));
  resid_dl_64.push_back(makeGraph("output/resid_dl_6arrays_e64.txt", "host a 64", kBlue, kOpenTriangleDown));
  //  resid_dl_8.push_back(makeGraph("output/resid_dl_6arrays_e8.txt", "host a 8", kBlue, kOpenTriangleDown));

  resid_dl_2097152.push_back(makeGraph("output_mic/resid_dl_shsa_e2097152.txt", "mic shsa 2097152", kOrange, kOpenTriangleUp));
  resid_dl_262144.push_back(makeGraph("output_mic/resid_dl_shsa_e262144.txt", "mic shsa 262144", kOrange, kOpenTriangleUp));
  resid_dl_32768.push_back(makeGraph("output_mic/resid_dl_shsa_e32768.txt", "mic shsa 32768", kOrange, kOpenTriangleUp));
  resid_dl_4096.push_back(makeGraph("output_mic/resid_dl_shsa_e4096.txt", "mic shsa 4096", kOrange, kOpenTriangleUp)); 
  resid_dl_512.push_back(makeGraph("output_mic/resid_dl_shsa_e512.txt", "mic shsa 512", kOrange, kOpenTriangleUp));
  resid_dl_64.push_back(makeGraph("output_mic/resid_dl_shsa_e64.txt", "mic shsa 64", kOrange, kOpenTriangleUp));
  //  resid_dl_8.push_back(makeGraph("output_mic/resid_dl_shsa_e8.txt", "mic shsa 8", kOrange, kOpenTriangleUp));
  
  resid_dl_2097152.push_back(makeGraph("output_mic/resid_dl_6arrays_e2097152.txt", "mic a 2098152", kCyan, kOpenSquare));
  resid_dl_262144.push_back(makeGraph("output_mic/resid_dl_6arrays_e262144.txt", "mic a 262144", kCyan, kOpenSquare));
  resid_dl_32768.push_back(makeGraph("output_mic/resid_dl_6arrays_e32768.txt", "mic a 32768", kCyan, kOpenSquare));
  resid_dl_4096.push_back(makeGraph("output_mic/resid_dl_6arrays_e4096.txt", "mic a 4096", kCyan, kOpenSquare)); 
  resid_dl_512.push_back(makeGraph("output_mic/resid_dl_6arrays_e512.txt", "mic a 512", kCyan, kOpenSquare));
  resid_dl_64.push_back(makeGraph("output_mic/resid_dl_6arrays_e64.txt", "mic a 64", kCyan, kOpenSquare));
  //  resid_dl_8.push_back(makeGraph("output_mic/resid_dl_6arrays_e8.txt", "mic a 8", kCyan, kOpenSquare));

  TCanvas *resid_dl_2098152_c = drawGraphs(resid_dl_2097152);
  TCanvas *resid_dl_262144_c = drawGraphs(resid_dl_262144);
  TCanvas *resid_dl_32768_c = drawGraphs(resid_dl_32768);
  TCanvas *resid_dl_4096_c = drawGraphs(resid_dl_4096);
  TCanvas *resid_dl_512_c = drawGraphs(resid_dl_512);
  TCanvas *resid_dl_64_c = drawGraphs(resid_dl_64);
  //  TCanvas *resid_dl_8_c = drawGraphs(resid_dl_8);


  //  resid_host_dl_array.push_back(makeGraph("output/resid_dl_6arrays_e2097152_v8_i262144.txt", "8", kRed, kCircle));
//   resid_host_dl_array.push_back(makeGraph("output/resid_dl_6arrays_e2097152_v16_i131072.txt", "16", kRed, kCircle));
  // resid_host_dl_array.push_back(makeGraph("output/resid_dl_6arrays_e2097152_v32_i65536.txt", "32", kRed, kCircle));
  // resid_host_dl_array.push_back(makeGraph("output/resid_dl_6arrays_e2097152_v64_i32768.txt", "64", kRed, kCircle));
  // resid_host_dl_array.push_back(makeGraph("output/resid_dl_6arrays_e2097152_v128_i16384.txt", "128", kRed, kCircle));
  // resid_host_dl_array.push_back(makeGraph("output/resid_dl_6arrays_e2097152_v256_i8192.txt", "256", kRed, kCircle));
  // resid_host_dl_array.push_back(makeGraph("output/resid_dl_6arrays_e2097152_v512_i4096.txt", "512", kRed, kCircle));
  // resid_host_dl_array.push_back(makeGraph("output/resid_dl_6arrays_e2097152_v1024_i2048.txt", "1024", kRed, kCircle));
  // resid_host_dl_array.push_back(makeGraph("output/resid_dl_6arrays_e2097152_v2048_i1024.txt", "2048", kRed, kCircle));
  // resid_host_dl_array.push_back(makeGraph("output/resid_dl_6arrays_e2097152_v4096_i512.txt", "4096", kRed, kCircle));
  //resid_host_dl_array.push_back(makeGraph("output/resid_dl_6arrays_e2097152_v8192_i256.txt", "8192", kRed, kCircle));
  //  resid_host_dl_array.push_back(makeGraph("output/resid_dl_6arrays_e2097152_v16384_i128.txt", "16384", kRed, kCircle));
  //resid_host_dl_array.push_back(makeGraph("output/resid_dl_6arrays_e2097152_v32768_i64.txt", "32768", kRed, kCircle));
  //resid_host_dl_array.push_back(makeGraph("output/resid_dl_6arrays_e2097152_v65536_i32.txt", "65536", kRed, kCircle));
  //  resid_host_dl_array.push_back(makeGraph("output/resid_dl_6arrays_e2097152_v131072_i16.txt", "131072", kRed, kCircle));
  //resid_host_dl_array.push_back(makeGraph("output/resid_dl_6arrays_e2097152_v262144_i8.txt", "262144", kRed, kCircle));
  //resid_host_dl_array.push_back(makeGraph("output/resid_dl_6arrays_e2097152_v524288_i4.txt", "524288", kRed, kCircle));
  //resid_host_dl_array.push_back(makeGraph("output/resid_dl_6arrays_e2097152_v1048576_i2.txt", "1048576", kRed, kCircle));
  //  resid_host_dl_array.push_back(makeGraph("output/resid_dl_6arrays_e2097152_v2097152_i1.txt", "2097152", kRed, kCircle));
  //  TCanvas *resid_host_dl_c = drawGraphs(resid_host_dl_array);










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
  temp->SetMarkerColor(color);
  temp->SetLineColor(color);
  temp->SetMarkerStyle(shape);
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
	graphs[0]->Draw("APC");
	l->AddEntry(graphs[0], graphs[0]->GetTitle(), "lp"); 
  }
  for(unsigned int i=1; i<graphs.size(); ++i){
	graphs[i]->Draw("same PC");
	l->AddEntry(graphs[i],graphs[i]->GetTitle(), "lp");
  }
  l->Draw();
  return c;
}



//}
