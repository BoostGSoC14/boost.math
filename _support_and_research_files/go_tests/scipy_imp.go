package main

// This implementation just executes python interpreter
// and calls functions. It is ineffective. Use cephes_scipy_imp
// instead

import (
  "os/exec"
  "bytes"
  "strconv"
  "strings"
)

func hyp1f1_python(a, b, z float64) (res string, err error){
  const python_cmd = "python"
  python_args := []string{"-c"}

  const import_cmd = "from mpmath import *;"

  a_str, b_str, z_str :=
    strconv.FormatFloat(a, 'e', 15, 64),
    strconv.FormatFloat(b, 'e', 15, 64),
    strconv.FormatFloat(z, 'e', 15, 64)

  full_cmd := strings.Join(
      []string{import_cmd, "print hyp1f1(", a_str, ",", b_str, ",", z_str, ")"},
      " ",
  )

  command_slice := append(python_args, full_cmd)

  out, err := exec.Command(python_cmd, command_slice...).Output()

  res = bytes.NewBuffer(out).String()
  return res, err
}
