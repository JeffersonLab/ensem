## Unittests for ensem library

import ensem, complex
import unittest
import random

suite "Test the ensem suite":
  const
    nbin  = 3
    Lt    = 4
    fluct = 2.0
    EPS   = 1.0e-8

  proc buildEnsem(avg: float64; max: float64): Ensemble_t =
    ## Build up a random number sequence about an average `avg`
    result = newEnsemble(RealType, nbin, Lt)
    for n in 0..nbin-1:
      # Some simple baseline
      var dest = newSeq[float64](Lt)
      for t in 0..Lt-1:
        dest[t] = avg + rand(max)

      result[n] = dest

  proc buildEnsem(avg: Complex64; max: float64): Ensemble_t =
    ## Build up a random number sequence about an average `avg`
    let re = buildEnsem(avg.re, max)
    let im = buildEnsem(avg.im, max)
    result = cmplx(re,im)

     
  test "Can we build a real ensemble, serialize it, and deserialize it":
    let x = buildEnsem(17.5, fluct)
    let xx = $x
    echo "A real ensemble:\n", xx
    let y = deserializeEnsemble(xx)
    require(x =~ y)

  test "Try calc-ing a real one":
    let x = buildEnsem(17.5, fluct)
    echo "Calc a real ensemble:"
    let vv = calc(x)
    echo vv
    let v: CalcDataType_t = vv.data[1]
    echo "Extract the avg/err from a time-slice: ", v.avg, " +- ", v.err
    let xx = $x
    let y = deserializeEnsemble(xx)
    require(x =~ y)
    # Warning - the tests here assume a certain order of calls to the RNG
    # If you move around these tests, they will break
    require(abs(v.avg.re - 18.68977840556422) < EPS)
    require(abs(v.err - 0.4563763365256296) < EPS)

  test "Try calc-ing a complex one":
    let x = buildEnsem(complex64(45.0f64, 80.5f64), fluct)
    echo "Calc a complex ensemble:\n", calc(x)
    let xx = $x
    let y = deserializeEnsemble(xx)
    require(x =~ y)

  test "Addition":
    let src1 = buildEnsem(14.5, fluct)
    let src2 = buildEnsem(42.0, fluct)
    let x = src1 + src2
    let xx = $x
    let y = deserializeEnsemble(xx)
    require(x =~ y)

  test "Multiplication":
    let src1 = buildEnsem(complex64(20.0, 5.6), fluct)
    let src2 = buildEnsem(42, fluct)
    let x = src1 * src2
    echo "calc the product of a complex and real ensemble:\n", calc(x)
    let xx = $x
    let y = deserializeEnsemble(xx)
    require(x =~ y)

  test "More complicated function":
    let src1 = buildEnsem(complex64(20.0, 5.6), fluct)
    let src2 = buildEnsem(4, fluct)
    let x = calc(real(src1 * exp(-src2)))
    echo "calc a more complicated expr:\n", x
    let v = x.data[0]
    # Warning - the tests here assume a certain order of calls to the RNG
    # If you move around these tests, they will break
    require(abs(v.avg.re - 0.1427419969424772) < EPS)
    require(abs(v.err - 0.04424750135381803) < EPS)
