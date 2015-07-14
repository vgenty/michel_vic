#include <vector>
#include <iostream>

int main() {

  std::vector<double> cluster_hits = {1,2,3,4,5,6,7,8,8,9,9};
  std::vector<double> michel_hits  = {1,10,12,13,14,2,3};
  
  std::cout << "cluster_hits: { ";
  for(const auto& c : cluster_hits)
    std::cout << c << ",";
  std::cout << "}";
  std::cout << std::endl;
  std::cout << "michel_hits: { ";
  for(const auto& m : michel_hits)
    std::cout << m << ",";
  std::cout << "}";
  std::cout << std::endl;
  
  std::vector<double>::iterator a1;
  int  w = 0;
  bool dup = true;
  int  ndupes = 0;
  
  while(1) {
    for(a1 = cluster_hits.begin(); a1 != cluster_hits.end(); ++a1) {
      for(const auto& m : michel_hits) {
	if(*a1 == m) {
	  w++;
	  break;  
	}
      }
      if(w > 0)
	break;
    }
     
    if(w > 0) {
      cluster_hits.erase(a1);
      w = 0;
      dup = true;
      ndupes++;
    }
    w = 0;
     
    if(!dup)
      break;
     
    dup = false;
  }

  std::cout << "encountered : " << ndupes << "\n";
  std::cout << "cluster_hits: { ";
  for(const auto& c : cluster_hits)
    std::cout << c << ",";
  std::cout << "}";
  std::cout << std::endl;
  
  
  std::vector<double>::iterator a2, a3;
  std::cout << "...checking if cluster_hits in cluster_hits..\n";
  w = 0;
  dup = false;
  int ndupes2 = 0;
  
  while(1) {
    for(a2 = cluster_hits.begin(); a2 != cluster_hits.end(); ++a2) {
      for(a3 = cluster_hits.begin(); a3 != cluster_hits.end(); ++a3) {
	if(*a2 == *a3) {
	  w++;
	}
	if(w > 1)
	  break;  
      }
      if(w > 1)
	break;
      w = 0;
    }
    if(w > 1) {
      cluster_hits.erase(a2);
      w = 0;
      dup = true;
      ndupes2++;
    }
    w = 0;
    
    if(!dup)
      break;
    
    dup = false;
  }
  
  std::cout << "encountered : " << ndupes + ndupes2<< "\n";
  std::cout << "cluster_hits: { ";
  for(const auto& c : cluster_hits)
    std::cout << c << ",";
  std::cout << "}";
  std::cout << std::endl;
  
  
}
