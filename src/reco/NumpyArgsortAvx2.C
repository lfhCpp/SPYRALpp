#include "x86simdsort-static-incl.h"
namespace rk4 {
void NumpyArgsortAvx2(double* values,size_t* indices,size_t n){
    x86simdsortStatic::argsort(values,indices,n,true);
}
}
