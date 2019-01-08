## Support for ensemble file format and arithmetic using jackknife/bootstrap propagation of errors.

import
  complex, math, system, streams, parseutils, strutils

type
  DataType_t* = enum
    RealType = 0,
    ComplexType = 1

  EnsemType_t* = enum
    EnsemJackknife = 10,
    EnsemBootstrap = 11

  CalcDataType_t* = tuple[avg: Complex64, err: float64]

  CalcType_t* = object
    data*:       seq[CalcDataType_t]
    typ*:        DataType_t      ## 0 if real and 1 if complex

  Ensemble_t* = object
    data:        seq[Complex64]    ## Data in either rescaled or raw format
    typ:         DataType_t             ## 0 if real and 1 if complex
    nbin:        int                    ## Number of bins, or data samples, configs, etc
    Lt:          int                    ## Number of time-slices of a prop
#    ens:         EnsemType_t           ## Type of fluctuations for the ensemble


proc norm2(x: Complex64): float =
  result = x.re * x.re + x.im * x.im

proc MAX(x: int; y: int): int =
  return if (x < y): y else: x

proc MIN(x: int; y: int): int =
  return if (x < y): x else: y


proc newCalcType*(typ: DataType_t; Lt: int): CalcType_t =
  ## Create a new calc object of length `Lt`
  result.typ  = typ
  result.data = newSeq[CalcDataType_t](Lt)


proc dataType*(src: CalcType_t): DataType_t = 
  ## Get the data-type
  result = src.typ


proc numElem*(src: Ensemble_t): int = 
  ## Get the time-extent
  result = src.Lt


proc numBins*(src: Ensemble_t): int = 
  ## Get the number of bins
  result = src.nbin


proc dataType*(src: Ensemble_t): DataType_t = 
  ## Get the data-type
  result = src.typ


proc promote_type(src1: DataType_t; src2: DataType_t): DataType_t =
  ## Promote type of data 
  result = RealType
  ## just to make compiler happy 
  if src1 == RealType and src2 == RealType:
    result = RealType
  elif src1 == RealType and src2 == ComplexType:
    result = ComplexType
  elif src1 == ComplexType and src2 == RealType:
    result = ComplexType
  elif src1 == ComplexType and src2 == ComplexType:
    result = ComplexType
  else:
    quit("some unknown types in promote_type")


proc rescale_factor(num: int; typ: EnsemType_t): float =
  ## Return the proper rescaling factor appropriate for the EnsemType of `typ`
  case typ:
    of EnsemJackknife:
      result = -(float(num - 1))
    of EnsemBootstrap:
      result = sqrt(float(num - 1))
#    else:
#      quit("unknown EnsemType= " & $typ)


proc rescaleFactor*(src: Ensemble_t): float =
  ## Return the proper rescaling factor appropriate for the ensemble
  result = rescale_factor(src.nbin, EnsemJackknife)


proc newEnsemble*(ens: EnsemType_t; typ: DataType_t; nbin: int; Lt: int): Ensemble_t =
  ## Create a new ensemble of some given number of bins and Lt 
  if nbin <= 0 or Lt <= 0:
    quit("invalid input in newEnsemble")
  result.data = newSeq[Complex64](nbin*Lt)
  if result.data.len == 0:
    quit("malloc returned NULL in newEnsemble")
  result.typ  = typ
  result.nbin = nbin
  result.Lt   = Lt
#  result.ens  = ens


proc newEnsemble*(typ: DataType_t; nbin: int; Lt: int): Ensemble_t =
  ## Create a new ensemble of some given number of bins and Lt 
  result = newEnsemble(EnsemJackknife, typ, nbin, Lt)


proc newEnsemble*(src: Ensemble_t): Ensemble_t =
  ## Create a new ensemble using parameters from source 
#  result = newEnsemble(src.ens, src.typ, src.nbin, src.Lt)
  result = newEnsemble(EnsemJackknife, src.typ, src.nbin, src.Lt)


proc newEnsemble*(val: SomeNumber; nbin: int; Lt: int): Ensemble_t =
  ## Promote constant to ensemble 
  result = newEnsemble(RealType, nbin, Lt)
  let v = float(val)
  var n = 0
  while n < nbin:
    var k = 0
    while k < Lt:
      result.data[k + n * Lt].re = v
      result.data[k + n * Lt].im = 0.0
      inc(k)
    inc(n)


proc newEnsemble*(val: Complex; nbin: int; Lt: int): Ensemble_t =
  ## Promote constant to ensemble 
  result = newEnsemble(ComplexType, nbin, Lt)
  var n = 0
  while n < nbin:
    var k = 0
    while k < Lt:
      result.data[k + n * Lt] = val
      inc(k)
    inc(n)


proc newEnsemble*(val: seq[SomeNumber]; nbin: int): Ensemble_t =
  ## Promote a sequence of numbers to an ensemble 
  let Lt = val.len
  result = newEnsemble(RealType, nbin, Lt)
  var n = 0
  while n < nbin:
    var k = 0
    while k < Lt:
      result.data[k + n * Lt].re = float(val[k])
      result.data[k + n * Lt].im = 0.0
      inc(k)
    inc(n)


proc newEnsemble*(val: seq[Complex]; nbin: int): Ensemble_t =
  ## Promote a sequence of numbers to an ensemble 
  let Lt = val.len
  result = newEnsemble(ComplexType, nbin, Lt)
  var n = 0
  while n < nbin:
    var k = 0
    while k < Lt:
      result.data[k + n * Lt] = val[k]
      inc(k)
    inc(n)


proc `[]`*[I: Natural](src: Ensemble_t; i: I): seq[Complex] =
  ## Return a particular bin
  if (i < 0) or (i >= src.nbin):
    quit("Index out-of-range of the ensemble: i= " & $i)
  let Lt = src.Lt
  result = newSeq[Complex](Lt)
  var k = 0
  while k < Lt:
    result[k] = src.data[k + i * Lt]
    inc(k)


proc `[]=`*[I: Natural](src: var Ensemble_t; i: I; val: seq[SomeNumber]) =
  ## Return a particular bin with a `val`. The `src` is expected to be a RealType ensemble.
  if src.typ != RealType:
    quit("Expected real type ensemble if assigning to a seq[numbers]")
  if (i < 0) or (i >= src.nbin):
    quit("Index out-of-range of the ensemble: i= " & $i)
  if val.len != src.Lt:
    quit("Non-commensurate time-extents: val.len = " & $val.len)
  let Lt = src.Lt
  var k = 0
  while k < Lt:
    src.data[k + i * Lt].re = float(val[k])
    src.data[k + i * Lt].im = 0.0
    inc(k)


proc `[]=`*[I: Natural](src: var Ensemble_t; i: I; val: seq[Complex]) =
  ## Return a particular bin with a `val`. The `src` is expected to be a ComplexType ensemble.
  if src.typ != ComplexType:
    quit("Expected complex type ensemble if assigning to a seq[numbers]")
  if (i < 0) or (i >= src.nbin):
    quit("Index out-of-range of the ensemble: i= " & $i)
  if val.len != src.Lt:
    quit("Non-commensurate time-extents: val.len = " & $val.len)
  let Lt = src.Lt
  var k = 0
  while k < Lt:
    src.data[k + i * Lt] = val[k]
    inc(k)



proc `$`*(z: CalcType_t): string = 
  ## Returns a nice string representation for a calc-ed object
  result = ""
  let Lt = z.data.len
  case z.typ
  of RealType:
    var k = 0
    while k < Lt:
      if k > 0: result &= "\n"
      let avg = z.data[k].avg
      let err = z.data[k].err
      if avg.re != 0.0:
        let rat = err / avg.re
        result &= $k & "   " & $avg.re & " " & $err & "   " & $rat
      else:
        result &= $k & "   " & $avg.re & " " & $err
      inc(k)

  of ComplexType:
    var k = 0
    while k < Lt:
      if k > 0: result &= "\n"
      let avg = z.data[k].avg
      let err = z.data[k].err
      var re = formatEng(avg.re, precision=12, trim= true)
      var im = formatEng(avg.im, precision=12, trim= true)
      result &= $k & "   ( " & re & " , " & im & " )   " & $err
      inc(k)

#  else:
#    quit("something wrong with calc type= " & $z.typ)


proc `$`*(src: Ensemble_t): string =
  ## Print the ensemble into a string
  let typ = src.typ
  let num = src.nbin
  let Lt  = src.Lt

  result = $num & " " & $Lt & " " & $int(typ) & " 0 1"
  case typ
  of RealType:
    for n in 0..num-1:
      for k in 0..Lt-1:
        var re = formatEng(src.data[k + Lt * n].re, precision=12)
        result &= "\n" & $k & " " & re

  of ComplexType:
    for n in 0..num-1:
      for k in 0..Lt-1:   
        var re = formatEng(src.data[k + Lt * n].re, precision=12, trim=true)
        var im = formatEng(src.data[k + Lt * n].im, precision=12, trim=true)
        result &= "\n" & $k & " " & re & " " & im

#  else:
#    quit("something wrong with ensemble: type= " & $typ)



proc check_two_ensemble(src1, src2: Ensemble_t): int =
  ## Check if two ensemble have the same parameters 
  if src1.nbin != src2.nbin:
    result = -1
  if src1.Lt == src2.Lt:
    result = 0
  elif src1.Lt == 1 or src2.Lt > 1:
    result = 1
  elif src1.Lt > 1 or src2.Lt == 1:
    result = 2
  else:
    result = -1


proc `=~`*[T](x, y: Complex[T]): bool =
  result = abs(x.re-y.re) < 1e-6 and abs(x.im-y.im) < 1e-6


proc `=~`*(src1, src2: CalcDataType_t): bool =
  ## Compare two ensembles `src1` and `src2` approximately
  result = (src1.avg =~ src2.avg)


proc `=~`*(src1, src2: Ensemble_t): bool =
  ## Compare two ensembles `src1` and `src2` approximately
  result = true
  if check_two_ensemble(src1, src2) != 0: return false
  let num = src1.nbin * src1.Lt
  for j in 0..num-1:
    let foo: bool = (src1.data[j] =~ src2.data[j])
    result = result and foo


proc rescale_ensemble(src: Ensemble_t; factor: float): Ensemble_t =
  ## Return a new rescaled ensemble 
  result = newEnsemble(src)
  let num = src.nbin
  let Lt  = src.Lt
  var avg: Complex64

  for k in 0..Lt-1:
    avg = complex64(0.0, 0.0)
    for n in 0..num-1:
      avg += src.data[k + Lt * n]

    avg = avg / float(num)
    for n in 0..num-1:
      result.data[k + Lt * n] = avg + (src.data[k + Lt * n] - avg) * factor


proc rescaleEnsemUp(src: Ensemble_t): Ensemble_t =
  ## Return a new ensemble with fluctuations rescaled upwards
  result = rescale_ensemble(src, src.rescale_factor)
  

proc rescaleEnsemDown(src: Ensemble_t): Ensemble_t =
  ## Return a new ensemble with fluctuations rescaled downwards
  result = rescale_ensemble(src, 1.0 / src.rescale_factor)
  


proc `+`*(val: SomeNumber; src2: Ensemble_t): Ensemble_t =
  ## Add constant on an ensemble 
  result = newEnsemble(src2)
  let v  = float(val)
  let Lt = src2.Lt
  for k in 0..Lt-1:
    for n in 0..src2.nbin-1:
      result.data[k + Lt * n] = v + src2.data[k + Lt * n]


proc `+`*(val: Complex; src2: Ensemble_t): Ensemble_t =
  ## Add constant on an ensemble 
  result = newEnsemble(src2)
  let Lt = src2.Lt
  for k in 0..Lt-1:
    for n in 0..src2.nbin-1:
      result.data[k + Lt * n] = val + src2.data[k + Lt * n]


proc `+`*(src1: Ensemble_t; val: SomeNumber): Ensemble_t =
  ## Add constant on an ensemble 
  result = val + src1


proc `+`*(src1: Ensemble_t; val: Complex): Ensemble_t =
  ## Add constant on an ensemble 
  result = val + src1


proc `+`*(src1: Ensemble_t; src2: Ensemble_t): Ensemble_t =
  ## Add two ensembles 
  let num = src1.nbin
  let Lt1 = src1.Lt
  let Lt2 = src2.Lt

  case check_two_ensemble(src1, src2)
  of 0, 2:
    result = newEnsemble(src1)
  of 1:
    result = newEnsemble(src2)
  else:
    quit("ensembles not compatible")

  result.typ = promote_type(src1.typ, src2.typ)
  let Lt = MAX(Lt1, Lt2)
  for k in 0..Lt-1:
    let k1 = MIN(k, Lt1 - 1)
    let k2 = MIN(k, Lt2 - 1)
    for n in 0..num-1:
      result.data[k + Lt * n] = src1.data[k1 + Lt1 * n] + src2.data[k2 + Lt2 * n]


proc `-`*(src1: Ensemble_t; src2: Ensemble_t): Ensemble_t =
  ## Subtract two ensembles 
  let num = src1.nbin
  let Lt1 = src1.Lt
  let Lt2 = src2.Lt
  var
    Lt: int
    n: int
    k: int
    k1: int
    k2: int
  case check_two_ensemble(src1, src2)
  of 0, 2:
    result = newEnsemble(src1)
  of 1:
    result = newEnsemble(src2)
  else:
    quit("ensembles not compatible")

  result.typ = promote_type(src1.typ, src2.typ)
  Lt = MAX(Lt1, Lt2)
  k = 0
  while k < Lt:
    k1 = MIN(k, Lt1 - 1)
    k2 = MIN(k, Lt2 - 1)
    n = 0
    while n < num:
      result.data[k + Lt * n] = src1.data[k1 + Lt1 * n] - src2.data[k2 + Lt2 * n]
      inc(n)
    inc(k)


proc `*`*(val: SomeNumber; src2: Ensemble_t): Ensemble_t =
  ## Multiply a constant on an ensemble 
  result = newEnsemble(src2)
  let v  = float(val)
  let Lt = src2.Lt
  for k in 0..Lt-1:
    for n in 0..src2.nbin-1:
      result.data[k + Lt * n] = v * src2.data[k + Lt * n]


proc `*`*(val: Complex; src2: Ensemble_t): Ensemble_t =
  ## Multiply a constant on an ensemble 
  result = newEnsemble(src2)
  let Lt = src2.Lt
  for k in 0..Lt-1:
    for n in 0..src2.nbin-1:
      result.data[k + Lt * n] = val * src2.data[k + Lt * n]


proc `*`*(src1: Ensemble_t; val: SomeNumber): Ensemble_t =
  ## Multiply a constant on an ensemble 
  result = val * src1


proc `*`*(src1: Ensemble_t; val: Complex): Ensemble_t =
  ## Multiply a constant on an ensemble 
  result = val * src1


proc `*`*(src1a: Ensemble_t; src2a: Ensemble_t): Ensemble_t =
  ## Multiply two ensembles 
  # First rescale downwards the fluctuations
  let src1 = rescaleEnsemDown(src1a)
  let src2 = rescaleEnsemDown(src2a)

  var dest: Ensemble_t
  let num = src1.nbin
  let Lt1 = src1.Lt
  let Lt2 = src2.Lt
  var
    Lt: int
    n: int
    k: int
    k1: int
    k2: int
  case check_two_ensemble(src1, src2)
  of 0, 2:
    dest = newEnsemble(src1)
  of 1:
    dest = newEnsemble(src2)
  else:
    quit("ensembles not compatible")

  dest.typ = promote_type(src1.typ, src2.typ)
  Lt = MAX(Lt1, Lt2)
  k = 0
  while k < Lt:
    k1 = MIN(k, Lt1 - 1)
    k2 = MIN(k, Lt2 - 1)
    n = 0
    while n < num:
      dest.data[k + Lt * n] = src1.data[k1 + Lt1 * n] * src2.data[k2 + Lt2 * n]
      inc(n)
    inc(k)
  # Rescale fluctuations upwards for the output
  result = rescaleEnsemUp(dest)


proc `/`*(src1a: Ensemble_t; src2a: Ensemble_t): Ensemble_t =
  ## Divide two ensembles 
  # First rescale downwards the fluctuations
  let src1 = rescaleEnsemDown(src1a)
  let src2 = rescaleEnsemDown(src2a)

  var dest: Ensemble_t
  let num = src1.nbin
  let Lt1 = src1.Lt
  let Lt2 = src2.Lt

  case check_two_ensemble(src1, src2)
  of 0, 2:
    dest = newEnsemble(src1)
  of 1:
    dest = newEnsemble(src2)
  else:
    quit("add: ensembles not compatible")

  dest.typ = promote_type(src1.typ, src2.typ)
  let Lt = MAX(Lt1, Lt2)

  for k in 0..Lt-1:
    let k1 = MIN(k, Lt1 - 1)
    let k2 = MIN(k, Lt2 - 1)
    ## For efficiency, split into real and complex 
    case src2.typ
    of RealType:
      for n in 0..num-1:
        dest.data[k + Lt * n].re = src1.data[k1 + Lt1 * n].re /
            src2.data[k2 + Lt2 * n].re
        dest.data[k + Lt * n].im = src1.data[k1 + Lt1 * n].im /
            src2.data[k2 + Lt2 * n].re

    of ComplexType:
      for n in 0..num-1:
        let denom: float = 1.0 /
            (src2.data[k2 + Lt2 * n].re * src2.data[k2 + Lt2 * n].re +
            src2.data[k2 + Lt2 * n].im * src2.data[k2 + Lt2 * n].im)
        dest.data[k + Lt * n].re = (src1.data[k1 + Lt1 * n].re *
            src2.data[k2 + Lt2 * n].re +
            src1.data[k1 + Lt1 * n].im * src2.data[k2 + Lt2 * n].im) * denom
        dest.data[k + Lt * n].im = (src1.data[k1 + Lt1 * n].im *
            src2.data[k2 + Lt2 * n].re -
            src1.data[k1 + Lt1 * n].re * src2.data[k2 + Lt2 * n].im) * denom

#    else:
#      quit("something wrong with src2 type: type= " & $src2.typ)

  # Rescale fluctuations upwards for the output
  result = rescaleEnsemUp(dest)


proc `-`*(src: Ensemble_t): Ensemble_t =
  ## Negate ensemble 
  result = newEnsemble(src)
  let Lt = src.Lt
  for k in 0..Lt-1:
    for n in 0..src.nbin-1:
      result.data[k + Lt * n] = - src.data[k + Lt * n]


proc replicate*(src: Ensemble_t; num: int): Ensemble_t =
  ## Replicate a src `num` times
  let N  = src.nbin * num
  let Lt = src.Lt
  result = newEnsemble(src.dataType, N, Lt)
  var j = 0
  while j < N:
    var n = 0
    while n < src.nbin*Lt:
      result.data[j] = src.data[n]
      inc(n)
      inc(j)


proc real*(src: Ensemble_t): Ensemble_t =
  ## Real part ensemble 
  result = newEnsemble(src)
  let Lt = src.Lt
  result.typ = RealType
  for k in 0..Lt-1:
    for n in 0..src.nbin-1:
      result.data[k + Lt * n].re = src.data[k + Lt * n].re
      result.data[k + Lt * n].im = 0.0


proc imag*(src: Ensemble_t): Ensemble_t =
  ## Imag part ensemble 
  result = newEnsemble(src)
  let Lt = src.Lt
  result.typ = RealType
  for k in 0..Lt-1:
    for n in 0..src.nbin-1:
      result.data[k + Lt * n].re = src.data[k + Lt * n].im
      result.data[k + Lt * n].im = 0.0


proc conjugate*(src: Ensemble_t): Ensemble_t = 
  ## Conjugate the ensemble 
  result = newEnsemble(src)
  let Lt = src.Lt

  ## Decide based on the input type 
  case src.typ
  of RealType:
    result.typ = RealType
    for k in 0..Lt-1:
      for n in 0..src.nbin-1:
        result.data[k + Lt * n].re = src.data[k + Lt * n].re
        result.data[k + Lt * n].im = 0.0

  of ComplexType:
    result.typ = ComplexType
    for k in 0..Lt-1:
      for n in 0..src.nbin-1:
        result.data[k + Lt * n].re = src.data[k + Lt * n].re
        result.data[k + Lt * n].im = - src.data[k + Lt * n].im

#  else:
#    quit("something wrong with src type: type= " & $src.typ)


proc norm2*(srca: Ensemble_t): Ensemble_t =
  ## Norm2 the ensemble 
  let src = rescaleEnsemDown(srca)
  var dest = newEnsemble(src)
  let Lt   = src.Lt

  ## Decide based on the input type 
  case src.typ
  of RealType:
    dest.typ = RealType
    for k in 0..Lt-1:
      for n in 0..src.nbin-1:
        let re: float = src.data[k + Lt * n].re
        dest.data[k + Lt * n].re = re * re
        dest.data[k + Lt * n].im = 0.0

  of ComplexType:
    dest.typ = RealType
    for k in 0..Lt-1:
      for n in 0..src.nbin-1:
        let re: float = src.data[k + Lt * n].re
        let im: float = src.data[k + Lt * n].im
        dest.data[k + Lt * n].re = re * re + im * im
        dest.data[k + Lt * n].im = 0.0

#  else:
#    quit("something wrong with src type: type= " & $src.typ)

  # Rescale fluctuations upwards for the output
  result = rescaleEnsemUp(dest)


proc arctan2*(src1a: Ensemble_t; src2a: Ensemble_t): Ensemble_t =
  ## Call atan2 on two real ensembles 
  # First rescale downwards the fluctuations
  let src1 = rescaleEnsemDown(src1a)
  let src2 = rescaleEnsemDown(src2a)

  var dest: Ensemble_t
  let num = src1.nbin
  let Lt1 = src1.Lt
  let Lt2 = src2.Lt

  case check_two_ensemble(src1, src2)
  of 0, 2:
    dest = newEnsemble(src1)
  of 1:
    dest = newEnsemble(src2)
  else:
    quit("cmplx: ensembles not compatible")

  if not (src1.typ == RealType and src2.typ == RealType):
    quit("atan2 requires both ensembles real")

  dest.typ = RealType
  let Lt = MAX(Lt1, Lt2)
  for k in 0..Lt-1:
    let k1 = MIN(k, Lt1 - 1)
    let k2 = MIN(k, Lt2 - 1)
    for n in 0..num-1:
      result.data[k + Lt * n].re = arctan2(src1.data[k1 + Lt1 * n].re,
                                   src2.data[k2 + Lt2 * n].re)
      result.data[k + Lt * n].im = 0.0


proc cmplx*(src1: Ensemble_t; src2: Ensemble_t): Ensemble_t =
  ## Build complex from two real ensembles 
  let num = src1.nbin
  let Lt1 = src1.Lt
  let Lt2 = src2.Lt

  case check_two_ensemble(src1, src2)
  of 0, 2:
    result = newEnsemble(src1)
  of 1:
    result = newEnsemble(src2)
  else:
    quit("cmplx: ensembles not compatible")

  if not (src1.typ == RealType and src2.typ == RealType):
    quit("cmplx requires both ensembles real")

  result.typ = ComplexType
  let Lt = MAX(Lt1, Lt2)
  for k in 0..Lt-1:
    let k1 = MIN(k, Lt1 - 1)
    let k2 = MIN(k, Lt2 - 1)
    for n in 0..num-1:
      result.data[k + Lt * n].re = src1.data[k1 + Lt1 * n].re
      result.data[k + Lt * n].im = src2.data[k2 + Lt2 * n].re


proc timesI*(src: Ensemble_t): Ensemble_t =
  ## Multiply ensemble by I 
  result  = newEnsemble(src)
  let num = src.nbin
  let Lt  = src.Lt

  result.typ = ComplexType
  case src.typ
  of RealType:
    for k in 0..Lt-1:
      for n in 0..num-1:
        result.data[k + Lt * n].re = 0.0
        result.data[k + Lt * n].im = src.data[k + Lt * n].re

  of ComplexType:
    for k in 0..Lt-1:
      for n in 0..num-1:
        result.data[k + Lt * n].re = - src.data[k + Lt * n].im
        result.data[k + Lt * n].im = src.data[k + Lt * n].re

#  else:
#    quit("something wrong with ensemble: type= " & $src.typ)


proc deserializeEnsemble*(input: string): Ensemble_t =
  ## Parse the string `input` as an ensemble and return it
  var
    num, Lt, ttype, junk, ncol: int
    typ: DataType_t

  # Turn the input string into a stringstream, and read away
  var fs = newStringStream(input) 

  var line = fs.readLine()
  if line == "":
    quit("Some error reading header from input string")
  else:
    let ll = splitWhiteSpace(line)
    if ll.len != 5:
      quit("Error in header line: " & line)

    num   = parseInt(ll[0])
    Lt    = parseInt(ll[1])
    ttype = parseInt(ll[2])
    junk  = parseInt(ll[3])
    ncol  = parseInt(ll[4])

  typ = cast[DataType_t](ttype)
  if ncol != 1:
    quit("only support 1 column in ensemble deserializer")

  if typ != RealType and typ != ComplexType:
    quit("error in datatype for file")

  result = newEnsemble(typ, num, Lt)

  var t, k: int
  var n = 0
  while n < num:
    case typ
    of RealType:
      k = 0
      while k < Lt: 
        if not fs.readLine(line):
          quit("error deserializing input string")
        #   echo "parse line= ", line
        let ll = splitWhiteSpace(line)
        if ll.len != 2:
          quit("Input needs space split values, got: " & line)

        t = parseInt(ll[0])
        if k != t:
          quit("error reading time slice data from input string: line= " & line)
          
        var dr: float
        discard parseFloat(ll[1], dr)
        result.data[k + Lt * n].re = dr
        result.data[k + Lt * n].im = 0.0
        inc(k)

    of ComplexType:
      k = 0
      while k < Lt: 
        if not fs.readLine(line):
          quit("error parsing input string = " & line)
        #   echo "parse line= ", line
        let ll = splitWhiteSpace(line)
        if ll.len != 3:
          quit("Input needs space split values, got: " & line)

        t = parseInt(ll[0])
        if k != t:
          quit("error parsing time slice data from line = " & line)
          
        var dr, di: float
        discard parseFloat(ll[1], dr)
        discard parseFloat(ll[2], di)
        result.data[k + Lt * n].re = dr
        result.data[k + Lt * n].im = di
        inc(k)
#    else:
#      quit("something wrong with ensemble: type= " & $typ)
    inc(n)

  fs.close()


proc readEnsemble*(filename: string): Ensemble_t =
  ## Read ensemble from `filename`
  # Slurp in the entire contents of the file
  result = deserializeEnsemble(readFile(filename))


proc writeEnsemble*(src: Ensemble_t; name: string) =
  ## Write ensemble 
  writeFile(name, $src)


proc apply_func_ensemble(funcptr: proc (x: float): float; srca: Ensemble_t): Ensemble_t =
  ## Apply function to ensemble 
  let src = rescaleEnsemDown(srca)
  var dest = newEnsemble(src)
  let num = src.nbin * src.Lt
  if src.typ != RealType:
    quit("only func(real) not supported")
  var j = 0
  while j < num:
    dest.data[j].re = funcptr(src.data[j].re)
    dest.data[j].im = 0.0
    inc(j)
  result = rescaleEnsemUp(dest)


proc log*(src: Ensemble_t): Ensemble_t =
  ## exp an ensemble
  result = apply_func_ensemble(math.ln, src)


proc log10*(src: Ensemble_t): Ensemble_t =
  ## exp an ensemble
  result = apply_func_ensemble(math.log10, src)


proc exp*(src: Ensemble_t): Ensemble_t =
  ## exp an ensemble
  result = apply_func_ensemble(math.exp, src)


proc sqrt*(src: Ensemble_t): Ensemble_t =
  ## sqrt an ensemble
  result = apply_func_ensemble(math.sqrt, src)


proc cos*(src: Ensemble_t): Ensemble_t =
  ## sqrt an ensemble
  result = apply_func_ensemble(math.cos, src)


proc sin*(src: Ensemble_t): Ensemble_t =
  ## sqrt an ensemble
  result = apply_func_ensemble(math.sin, src)


proc tan*(src: Ensemble_t): Ensemble_t =
  ## sqrt an ensemble
  result = apply_func_ensemble(math.tan, src)


proc acos*(src: Ensemble_t): Ensemble_t =
  ## sqrt an ensemble
  result = apply_func_ensemble(math.arccos, src)


proc asin*(src: Ensemble_t): Ensemble_t =
  ## sqrt an ensemble
  result = apply_func_ensemble(math.arcsin, src)


proc atan*(src: Ensemble_t): Ensemble_t =
  ## sqrt an ensemble
  result = apply_func_ensemble(math.arctan, src)


proc cosh*(src: Ensemble_t): Ensemble_t =
  ## sqrt an ensemble
  result = apply_func_ensemble(math.cosh, src)


proc sinh*(src: Ensemble_t): Ensemble_t =
  ## sqrt an ensemble
  result = apply_func_ensemble(math.sinh, src)


proc tanh*(src: Ensemble_t): Ensemble_t =
  ## sqrt an ensemble
  result = apply_func_ensemble(math.tanh, src)


proc pow*(srca: Ensemble_t; val: SomeNumber): Ensemble_t =
  ## Apply pow(src,const) to ensemble 
  let src = rescaleEnsemDown(srca)
  var dest = newEnsemble(src)
  let v  = float(val)
  if src.typ != RealType:
    quit("only func(real) not supported")

  for j in 0..src.data.len-1:
    dest.data[j].re = pow(src.data[j].re, v)
    dest.data[j].im = 0.0
  result = rescaleEnsemUp(dest)
 

proc calc_real_ensemble(src: Ensemble_t): CalcType_t =
  ## Calculate mean, err and err/mean for data 
  result  = newCalcType(RealType, src.Lt)
  let num = src.nbin
  let Lt  = src.Lt
  var avg: float
  var err: float

  for k in 0..Lt-1:
    avg = 0.0
    err = 0.0
    for n in 0..num-1:
      avg += src.data[k + Lt * n].re

    avg /= float(num)
    for n in 0..num-1:
      let diff = src.data[k + Lt * n].re - avg
      err += diff * diff

    err = sqrt(err / float((num - 1) * num))
    result.data[k] = (complex64(avg,0.0), err)

proc calc_complex_ensemble(src: Ensemble_t): CalcType_t =
  ## Calculate mean, err and err/mean for data 
  result = newCalcType(ComplexType, src.Lt)
  let num = src.nbin
  let Lt  = src.Lt
  var avg: Complex64
  var err: float

  for k in 0..Lt-1:
    avg = complex64(0.0, 0.0)
    err = 0.0
    for n in 0..num-1:
      avg += src.data[k + Lt * n]

    avg = avg / float(num)
    for n in 0..num-1:
      let diff = src.data[k + Lt * n] - avg
      err += norm2(diff)

    err = sqrt(err / float((num - 1) * num));
    result.data[k] = (avg, err)


proc calc*(src: Ensemble_t): CalcType_t =
  ## Calculate mean, err and err/mean for data 
  case src.typ
  of RealType:
    result = calc_real_ensemble(src)
  of ComplexType:
    result = calc_complex_ensemble(src)
#  else:
#    quit("something wrong with ensemble: type= " & $src.typ)


proc shift*(src: Ensemble_t; sh: int): Ensemble_t =
  ## Shift an ensemble in some direction dropping bits from the end 
  result  = newEnsemble(src)
  let num = src.nbin
  let Lt  = src.Lt
  var
    n: int
    k: int
    kk: int
  if abs(sh) > Lt:
    quit("shift: do not allow the shift greater than the Lt")

  n = 0
  while n < num:
    k = 0
    while k < Lt:
      kk = (k + Lt + sh) mod Lt
      result.data[k + Lt * n].re = src.data[kk + Lt * n].re
      result.data[k + Lt * n].im = src.data[kk + Lt * n].im
      inc(k)
    ## Clean out the ends 
    if sh > 0:
      k = Lt - sh
      while k < Lt:
        result.data[k + Lt * n].re = 0.0
        result.data[k + Lt * n].im = 0.0
        inc(k)
    elif sh < 0:
      k = 0
      while k < - sh:
        result.data[k + Lt * n].re = 0.0
        result.data[k + Lt * n].im = 0.0
        inc(k)
    inc(n)


proc cshift*(src: Ensemble_t; sh: int): Ensemble_t =
  ## Periodic (circular) shift an ensemble in some direction 
  result  = newEnsemble(src)
  let num = src.nbin
  let Lt  = src.Lt
  var
    n: int
    k: int
    kk: int
  n = 0
  while n < num:
    k = 0
    while k < Lt:
      kk = (k + Lt + sh) mod Lt
      result.data[k + Lt * n].re = src.data[kk + Lt * n].re
      result.data[k + Lt * n].im = src.data[kk + Lt * n].im
      inc(k)
    inc(n)


proc extract*(src: Ensemble_t; elem_i, elem_f: int): Ensemble_t =
  ## Extract a range of time slices from an ensemble 
  result.typ = src.typ
  result.typ = src.typ
  let num = src.nbin
  let Lt  = src.Lt
  var dLt = elem_f - elem_i
  var
    k: int
    n: int
  if elem_i < 0 or elem_i >= Lt:
    quit("index element out of bounds of ensemble: " & $elem_i)

  if elem_f < 0 or elem_f >= Lt:
    quit("index element out of bounds of ensemble: " & $elem_f)

  if elem_f < elem_i:
    quit("index element out of order in ensemble: " & $elem_i & " " & $elem_f)

  dLt = elem_f - elem_i + 1
  result = newEnsemble(src.typ, num, dLt)
  n = 0
  while n < num:
    k = 0
    while k < dLt:
      result.data[k + dLt * n].re = src.data[k + elem_i + Lt * n].re
      result.data[k + dLt * n].im = src.data[k + elem_i + Lt * n].im
      inc(k)
    inc(n)


proc extract*(src: Ensemble_t; elem_i: int): Ensemble_t =
  ## Extract a single time slice from an ensemble 
  result = extract(src, elem_i, elem_i)


proc merge*(src1: Ensemble_t; src2: Ensemble_t): Ensemble_t =
  ## Merge two ensembles in the observable time. Number of bins must be the same.
  if src1.nbin != src2.nbin:
    quit("Ensembles not compatible for merging")

  let Lt = src1.Lt + src2.Lt
  let num = src1.nbin
  let typ = promote_type(src1.typ, src2.typ)
  result = newEnsemble(typ, num, Lt)
  var n = 0
  while n < num:
    var k = 0
    while k < src1.Lt:
      result.data[k + n * Lt].re = src1.data[k + n * src1.Lt].re
      result.data[k + n * Lt].im = src1.data[k + n * src1.Lt].im
      inc(k)
    k = 0
    while k < src2.Lt:
      result.data[k + src1.Lt + n * Lt].re = src2.data[k + n * src2.Lt].re
      result.data[k + src1.Lt + n * Lt].im = src2.data[k + n * src2.Lt].im
      inc(k)
    inc(n)


proc concatMeas*(src1: Ensemble_t; src2: Ensemble_t): Ensemble_t =
  ## Append two ensembles 
  if (src1.Lt != src2.Lt) or (src1.typ != src2.typ):
    quit("Ensembles not compatible for appending")

  let Lt  = src1.Lt
  let num = src1.nbin + src2.nbin
  result  = newEnsemble(src1.typ, num, Lt)
  var n = 0
  while n < num:
    var j = 0
    while j < src1.nbin:
      var k = 0
      while k < Lt:
        result.data[k + n * Lt] = src1.data[k + j * Lt]
        inc(k)
      inc(j)
    j = 0
    while j < src2.nbin:
      var k = 0
      while k < Lt:
        result.data[k + n * Lt] = src2.data[k + j * Lt]
        inc(k)
      inc(j)
    inc(n)


proc appendMeas*(src: var Ensemble_t; src2: Ensemble_t) =
  ## Append an ensemble `src2` onto a `src` ensemble
  if (src.Lt != src2.Lt) or (src.typ != src2.typ):
    quit("Ensembles not compatible for appending")

  let Lt  = src.Lt
  let num = src.nbin
  src.nbin += src2.nbin
  src.data.setLen(src.nbin*Lt)
  var j = num
  var n = 0
  while n < src2.nbin:
    var k = 0
    while k < Lt:
      src.data[k + j * Lt] = src2.data[k + n * Lt]
      inc(k)
    inc(n)
    inc(j)


proc dropMeasEnd*(src: Ensemble_t; ndrop: Natural): Ensemble_t =
  ## Drop the last `ndrop` measurements from the `src` ensemble
  if src.nbin <= ndrop:
    quit("Dropping too many measurements: ndrop= " & $ndrop)
  let Lt  = src.Lt
  let num = src.nbin - ndrop
  result = newEnsemble(src.typ, num, Lt)
  var n = 0
  while n < num:
    var k = 0
    while k < Lt:
      result.data[k + n * Lt] = src.data[k + n * Lt]
      inc(k)
    inc(n)

