// #ifndef MALLOC_DEFS_H
// #define MALLOC_DEFS_H

// #include <iostream>
// #include <fstream>
// #include <sstream>
// #include <string>
// #include <cmath>
// #include <complex>
// #include <ginac/ginac.h>
// #include <cln/cln.h>
#include "mpc.h"
#include "global_vars.h"

using namespace std;


template<typename type>
extern void malloc_rk2_tens(type **&tens, int dim1, int dim2) {
  tens = new type*[dim1];
  for (int i=0; i<dim1; i++) {
    tens[i] = new type[dim2];
  }
}

template<typename type>
extern void mall_rk2_tens(type **&tens, int dim1, int dim2) {
  tens = (type**) malloc(dim1*sizeof(type*));
  for (int i=0; i<dim1; i++) {
    tens[i] = (type*) malloc(dim2*sizeof(type));;
  }
}

template<typename type>
void del_rk2_tens(type **&tens, int dim1) {
  for (int i=0; i<dim1; i++) {
    delete[] tens[i];
  }
  delete[] tens;
}

template<typename type>
void malloc_rk3_tens(type ***&tens, int dim1, int dim2, int dim3) {
  tens = new type**[dim1];
  for (int i=0; i<dim1; i++) {
    tens[i] = new type*[dim2];
    for (int j=0; j<dim2; j++) {
      tens[i][j] = new type[dim3];
    }
  }
}


template<typename type>
void malloc_rk4_tens(type ****&tens, int dim1, int dim2, int dim3, int dim4) {
  tens = new type***[dim1];
  for (int i=0; i<dim1; i++) {
    tens[i] = new type**[dim2];
    for (int j=0; j<dim2; j++) {
      tens[i][j] = new type*[dim3];
      for (int k=0; k<dim3; k++) {
        tens[i][j][k] = new type[dim4];
      }
    }
  }
}


template<typename type>
void malloc_rk5_tens(type ****&tens, int dim1, int dim2, int dim3, int dim4, int dim5) {
  tens = new type****[dim1];
  for (int i=0; i<dim1; i++) {
    tens[i] = new type***[dim2];
    for (int j=0; j<dim2; j++) {
      tens[i][j] = new type**[dim3];
      for (int k=0; k<dim3; k++) {
        tens[i][j][k] = new type*[dim4];
        for (int h=0; h<dim4; h++) {
          tens[i][j][k][h] = new type[dim5];
        }
      }
    }
  }
}


template<typename type>
void del_rk3_tens(type ***&tens, int dim1, int dim2) {
  for (int i=0; i<dim1; i++) {
    for (int j=0; j<dim2; j++) {
      delete[] tens[i][j];
    }
    delete[] tens[i];
  }
  delete[] tens;
}


template<typename type>
void del_rk4_tens(type ***&tens, int dim1, int dim2, int dim3) {
  for (int i=0; i<dim1; i++) {
    for (int j=0; j<dim2; j++) {
      for (int k=0; k<dim3; k++) {
        delete[] tens[i][j][k];
      }
      delete[] tens[i][j];
    }
    delete[] tens[i];
  }
  delete[] tens;
}


template<typename type>
void del_rk5_tens(type ***&tens, int dim1, int dim2, int dim3, int dim4) {
  for (int i=0; i<dim1; i++) {
    for (int j=0; j<dim2; j++) {
      for (int k=0; k<dim3; k++) {
        for (int h=0; h<dim4; h++) {
          delete[] tens[i][j][k][h];
        }
        delete[] tens[i][j][k];
      }
      delete[] tens[i][j];
    }
    delete[] tens[i];
  }
  delete[] tens;
}

// #endif