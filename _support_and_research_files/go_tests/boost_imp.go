package main

// #include "boost/instantiated.hpp"
// #cgo LDFLAGS: -L./boost -lboost_hyp1f1
import "C"

func hyp1f1_interface(a, b, z float64) float64{
  return float64(C.hypergeometric_1f1_interface_d(
        C.double(a),
        C.double(b),
        C.double(z),
        ),
  )
}

func hyp1f1_series(a, b, z C.double) float64{
  return float64(C.hypergeometric_1f1_series_d(
        C.double(a),
        C.double(b),
        C.double(z),
        ),
  )
}

func hyp1f1_asym(a, b, z C.double) float64{
  return float64(C.hypergeometric_1f1_asym_d(
        C.double(a),
        C.double(b),
        C.double(z),
        ),
  )
}

func hyp1f1_rational(a, b, z C.double) float64{
  return float64(C.hypergeometric_1f1_rational_d(
        C.double(a),
        C.double(b),
        C.double(z),
        ),
  )
}

func hyp1f1_13_3_7(a, b, z C.double) float64{
  return float64(C.hypergeometric_1f1_13_3_7_d(
        C.double(a),
        C.double(b),
        C.double(z),
        ),
  )
}
