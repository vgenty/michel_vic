#include <iostream>
#include <vector>


template<typename T>

std::vector<std::vector<T> > get_windows(const std::vector<T>& the_thing, int window_size)
{
  
  std::vector<std::vector<T> > data;
  std::vector<T> inner;
  auto w = window_size + 2;
  w = (unsigned int)((w - 1)/2);
  auto num = the_thing.size();
  
  for(int i = 1; i <= num; ++i) {
    if(i < w) {
      for(int j = 0; j < 2 * (i%w) - 1; ++j)
	inner.push_back(the_thing[j]);
    }else if (i > num - w + 1){
      for(int j = num - 2*((num - i)%w)-1 ; j < num; ++j)
	inner.push_back(the_thing[j]);
    }else{
      for(int j = i - w; j < i + w - 1; ++j)
	inner.push_back(the_thing[j]);
    }
    data.push_back(inner);
    inner.clear();
  }

  return data;
    
}

int main() {
  
  
  std::vector<size_t> aa = {1,2,3,4,5,6,7,8,9,10,11,12,13};
  std::cout << "aa : {";

  for(const auto& a : aa)
    std::cout << a << ",";
  std::cout << "}\n\n";
  
  auto bull = get_windows(aa,5);

  std::cout << "bull....\n {";
  for(const auto& b1: bull) {
    std::cout << "{";
    for(const auto& b2: b1) {
      std::cout << b2 << ",";
    }
    std::cout << "},";
  }
  std::cout << "}\n";

  
}
