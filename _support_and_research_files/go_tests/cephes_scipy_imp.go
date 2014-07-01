package main

// #cgo CFLAGS: -I/usr/include/python2.7
// #cgo LDFLAGS: -Lcephes_scipy/ -lcephes_scipy -lpython2.7
// #include "cephes_scipy/specfun_wrappers.h"
import "C"

func hyp1f1_cephes_scipy(a, b, z float64) float64 {
  return float64(C.hyp1f1_wrap(
        C.double(a),
        C.double(b),
        C.double(z),
        ),
  )
}
