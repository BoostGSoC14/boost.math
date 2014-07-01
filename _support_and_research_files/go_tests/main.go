package main

import (
  "fmt"
  "math/rand"
  "time"
)

func main() {
  r := rand.New(rand.NewSource(time.Now().UnixNano()))

  var a, b, c float64 =
    r.NormFloat64(),
    r.NormFloat64(),
    r.NormFloat64()

  fmt.Println(hyp1f1_interface(a, b, c))
  fmt.Println(hyp1f1_cephes_scipy(a, b, c))
}
