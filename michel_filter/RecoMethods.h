//singleton design pattern copied from who gives a fuck

#ifndef RECOMETHODS_H
#define RECOMETHODS_H
#include <iostream>

class RecoMethods
{
public:
  static RecoMethods& getInstance()
  {
    static RecoMethods    instance; // Guaranteed to be destroyed.
    // Instantiated on first use.
    return instance;
  }
  
  void aho() { std::cout << "bakakaka\n"; }
  
private:
  RecoMethods() {};                   // Constructor? (the {} brackets) are needed here.
  
  // // C++ 03
  // // ========
  // // Dont forget to declare these two. You want to make sure they
  // // are unacceptable otherwise you may accidentally get copies of
  // // your singleton appearing.
  // RecoMethods(RecoMethods const&);              // Don't Implement
  // void operator=(RecoMethods const&); // Don't implement
  
  // C++ 11
  // =======
  // We can use the better technique of deleting the methods
  // we don't want.

  RecoMethods(RecoMethods const&)     = delete;
  void operator=(RecoMethods const&)  = delete;
  
};

#endif

// /**
//  * \file RecoMethods.h
//  *
//  * \ingroup michel_filter
//  * 
//  * \brief Class def header for a class RecoMethods
//  *
//  * @author vgenty
//  */

// /** \addtogroup michel_filter

//     @{*/
// #ifndef RECOMETHODS_H
// #define RECOMETHODS_H

// #include <iostream>

// /**
//    \class RecoMethods
//    User defined class RecoMethods ... these comments are used to generate
//    doxygen documentation!
//  */
// class RecoMethods{

// public:

//   /// Default constructor
//   RecoMethods(){}

//   /// Default destructor
//   ~RecoMethods(){}

// };

// #endif
// /** @} */ // end of doxygen group 

