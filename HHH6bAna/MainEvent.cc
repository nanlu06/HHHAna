#define MainEvent_cxx

#include "MainEvent.h"
#include <TH2.h>
#include <TStyle.h>
#include<iostream>

bool MainEvent::mf_descend(int i,int j) { 
    return (i>j); 
}

double MainEvent::DeltaPhi(double phi1, double phi2) {
  double result = phi1 - phi2;
  while (result > M_PI)    result -= 2 * M_PI;
  while (result <= -M_PI)  result += 2 * M_PI;
  return result;
}

double MainEvent::DeltaR(double eta1, double phi1, double eta2, double phi2) {
  double deta = eta1 - eta2;
  double dphi = DeltaPhi(phi1, phi2);
  return std::sqrt(deta*deta + dphi*dphi);
}

//https://www.geeksforgeeks.org/keep-track-of-previous-indexes-after-sorting-a-vector-in-c-stl/
vector<pair<double, int>> MainEvent::sortArr(vector<double> arr, int n)
{
  
    // Vector to store element
    // with respective present index
    vector<pair<double, int> > vp;
  
    // Inserting element in pair vector
    // to keep track of previous indexes
    for (int i = 0; i < n; ++i) {
        vp.push_back(make_pair(arr[i], i));
    }
  
    // Sorting pair vector
    sort(vp.begin(), vp.end());
  
    // Displaying sorted element
    // with previous indexes
    // corresponding to each element
    /*
    cout << "Element\t"
         << "index" << endl;
    for (int i = 0; i < vp.size(); i++) {
        cout << vp[i].first << "\t"
             << vp[i].second << endl;
    }
    */
    return vp;
}
