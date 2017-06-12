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

  CalcDataType_t* = tuple[avg: Complex, err: float]

  CalcType_t* = object
    data*:       seq[CalcDataType_t]
    typ*:        DataType_t      ## 0 if real and 1 if complex

  Ensemble_t* = object
    data:        seq[Complex]    ## Data in either rescaled or raw format
    typ:         DataType_t      ## 0 if real and 1 if complex
    nbin:        int             ## Number of bins, or data samples, configs, etc
    Lt:          int             ## Number of time-slices of a prop
#    ens:         EnsemType_t    ## Type of fluctuations for the ensemble


proc norm2(x: Complex): float =
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
    quit("some unknown types in concatenate")


proc rescale_factor(num: int; typ: EnsemType_t): float =
  ## Return the proper rescaling factor appropriate for the EnsemType of `typ`
  case typ:
    of EnsemJackknife:
      result = -(float(num - 1))
    of EnsemBootstrap:
      result = sqrt(float(num - 1))
    else:
      quit("unknown EnsemType= " & $typ)


proc rescaleFactor*(src: Ensemble_t): float =
  ## Return the proper rescaling factor appropriate for the ensemble
  result = rescale_factor(src.nbin, EnsemJackknife)


proc newEnsemble*(ens: EnsemType_t; typ: DataType_t; nbin: int; Lt: int): Ensemble_t =
  ## Create a new ensemble of some given number of bins and Lt 
  if nbin <= 0 or Lt <= 0:
    quit("invalid input in newEnsemble")
  result.data = newSeq[Complex](nbin*Lt)
  if result.data == nil:
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
  let v = float64(val)
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
      result.data[k + n * Lt].re = float64(val[k])
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


proc `$`*(z: CalcType_t): string = 
  ## Returns a nice string representation for a calc-ed object
  result = ""
  let Lt = z.data.len
  if z.typ == RealType:
    var k = 0
    while k < Lt:
      let avg = z.data[k].avg
      let err = z.data[k].err
      if avg.re != 0.0:
        let rat = err / avg.re
        result &= $k & "   " & $avg.re & " " & $err & "   " & $rat & "\n"
      else:
        result &= $k & "   " & $avg.re & " " & $err & "   " & "\n"
      inc(k)

  if z.typ == ComplexType:
    var k = 0
    while k < Lt:
      let avg = z.data[k].avg
      let err = z.data[k].err
      var re = formatEng(avg.re, precision=12, trim= false)
      var im = formatEng(avg.im, precision=12, trim= false)
      result &= $k & "   ( " & re & " , " & im & " )   " & $err & "\n"
      inc(k)

  else:
    quit("something wrong with calc type= " & $z.typ)


proc `$`*(src: Ensemble_t): string =
  ## Print the ensemble into a string
  let typ = src.typ
  let num = src.nbin
  let Lt  = src.Lt

  result = $num & " " & $Lt & " " & $int(typ) & " 0 1" & "\n"
  case typ
  of RealType:
    for n in 0..num-1:
      for k in 0..Lt-1:
        var re = formatEng(src.data[k + Lt * n].re, precision=12)
        result &= $k & " " & re & "\n"

  of ComplexType:
    for n in 0..num-1:
      for k in 0..Lt-1:   
        var re = formatEng(src.data[k + Lt * n].re, precision=12, trim=false)
        var im = formatEng(src.data[k + Lt * n].im, precision=12, trim=false)
        result &= $k & " " & re & " " & im & "\n"

  else:
    quit("something wrong with ensemble: type= " & $typ)



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


proc `=~`*(src1, src2: Ensemble_t): bool =
  ## Compare two ensembles `src1` and `src2` approximately
  result = true
  if check_two_ensemble(src1, src2) != 0: return false
  let num = src1.nbin * src1.Lt
  for j in 0..num-1:
    let foo: bool = (src1.data[j] =~ src1.data[j])
    result = result and foo


proc rescale_ensemble(src: Ensemble_t; factor: float): Ensemble_t =
  ## Return a new rescaled ensemble 
  result = newEnsemble(src)
  let num = src.nbin
  let Lt  = src.Lt
  var avg: Complex

  for k in 0..Lt-1:
    avg = (0.0, 0.0)
    for n in 0..num-1:
      avg += src.data[k + Lt * n]

    avg /= float(num)
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
  let v  = float64(val)
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
  let v  = float64(val)
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

    else:
      quit("something wrong with src2 type: type= " & $src2.typ)

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

  else:
    quit("something wrong with src type: type= " & $src.typ)


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

  else:
    quit("something wrong with src type: type= " & $src.typ)

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

  else:
    quit("something wrong with ensemble: type= " & $src.typ)


proc readEnsemble*(name: string): Ensemble_t =
  ## Read ensemble 
  var
    num, Lt, ttype, junk, ncol: int
    typ: DataType_t

  # Slurp in the entire contents of the file
  var fs = newStringStream(readFile(name))
  if isNil(fs):
    quit("Error opening file = " & name)

  var line = fs.readLine()
  if line == "":
    quit("Some error reading header in file = " & name)
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
    quit("only support 1 column in " & name)

  if typ != RealType and typ != ComplexType:
    quit("error in type for file " & name)

  result = newEnsemble(typ, num, Lt)

  var t, k: int
  var n = 0
  while n < num:
    case typ
    of RealType:
      k = 0
      while k < Lt: 
        if not fs.readLine(line):
          quit("error reading file= " & name)
        #   echo "parse line= ", line
        let ll = splitWhiteSpace(line)
        if ll.len != 2:
          quit("Input needs space split values, got: " & line)

        t = parseInt(ll[0])
        if k != t:
          quit("error reading time slice data from " & name)
          
        var dr: float
        discard parseFloat(ll[1], dr)
        result.data[k + Lt * n].re = dr
        result.data[k + Lt * n].im = 0.0
        inc(k)

    of ComplexType:
      k = 0
      while k < Lt: 
        if not fs.readLine(line):
          quit("error reading file= " & name)
        #   echo "parse line= ", line
        let ll = splitWhiteSpace(line)
        if ll.len != 3:
          quit("Input needs space split values, got: " & line)

        t = parseInt(ll[0])
        if k != t:
          quit("error reading time slice data from " & name)
          
        var dr, di: float
        discard parseFloat(ll[1], dr)
        discard parseFloat(ll[2], di)
        result.data[k + Lt * n].re = dr
        result.data[k + Lt * n].im = di
        inc(k)
    else:
      quit("something wrong with ensemble: type= " & $typ)
    inc(n)

  fs.close()


proc writeEnsemble*(src: Ensemble_t; name: string) =
  ## Write ensemble 
  writeFile(name, $src)



proc apply_func_ensemble(funcptr: proc (x: float): float; src: Ensemble_t): Ensemble_t =
  ## Apply function to ensemble 
  result = newEnsemble(src)
  let num = src.nbin * src.Lt
  if src.typ != RealType:
    quit("only func(real) not supported")
  var j = 0
  while j < num:
    result.data[j].re = funcptr(src.data[j].re)
    result.data[j].im = 0.0
    inc(j)


proc exp*(src: Ensemble_t): Ensemble_t =
  ## exp an ensemble
  result = apply_func_ensemble(exp, src)


proc sqrt*(src: Ensemble_t): Ensemble_t =
  ## sqrt an ensemble
  result = apply_func_ensemble(math.sqrt, src)


proc pow*(src: Ensemble_t; val: SomeNumber): Ensemble_t =
  ## Apply pow(src,const) to ensemble 
  result = newEnsemble(src)
  let v  = float64(val)
  if src.typ != RealType:
    quit("only func(real) not supported")

  for j in 0..src.data.len-1:
    result.data[j].re = pow(src.data[j].re, val)
    result.data[j].im = 0.0
 

proc calc_real_ensemble(src: Ensemble_t): CalcType_t =
  ## Calculate mean, err and err/mean for data 
  result  = newCalcType(RealType, src.Lt)
  let num = src.nbin
  let Lt  = src.Lt

  for k in 0..Lt-1:
    var avg = 0.0
    var err = 0.0
    for n in 0..num-1:
      avg += src.data[k + Lt * n].re

    avg /= float(num)
    for n in 0..num-1:
      let diff = (src.data[k + Lt * n].re - avg)
      err += diff * diff

    err = sqrt(err / float((num - 1) * num))
    result.data[k] = ((avg,0.0), err)
#    if avg != 0.0: rat = err / avg
#    else: rat = 0.0
#    echo k, "   ", avg, " ", err, "   ", rat



proc calc_complex_ensemble(src: Ensemble_t): CalcType_t =
  ## Calculate mean, err and err/mean for data 
  result = newCalcType(ComplexType, src.Lt)
  let num = src.nbin
  let Lt  = src.Lt
  var avg: Complex
  var err: float

  for k in 0..Lt-1:
    avg = (0.0, 0.0)
    err = 0.0
    for n in 0..num-1:
      avg += src.data[k + Lt * n]

    avg /= float(num)
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
  else:
    quit("something wrong with ensemble: type= " & $src.typ)


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


proc concatenate*(src1: Ensemble_t; src2: Ensemble_t): Ensemble_t =
  ## Concatenate two ensembles 
  var
    k: int
    n: int
  var
    Lt: int
    num: int
  var typ: DataType_t
  if src1.nbin != src2.nbin:
    quit("Ensembles not compatible for concatenation")

  Lt = src1.Lt + src2.Lt
  num = src1.nbin
  typ = promote_type(src1.typ, src2.typ)
  result = newEnsemble(typ, num, Lt)
  n = 0
  while n < num:
    k = 0
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


#--------------------------------------------------------------------------
when isMainModule:
  let Lt   = 4
  let nbin = 3
  let src2 = newEnsemble(17.0, nbin, Lt)
  let src3 = newEnsemble((5.0,7.3), nbin, Lt)

  # Build a list in time
  var t = newSeq[float](Lt)
  for i in 0..Lt-1:
    t[i] = float(i)

  let src1 = newEnsemble(t, nbin)

#[
  echo "src1:"
  echo src1
  echo "src2:"
  echo src2
  echo "src3:"
  echo src3
]#

  var fred = src3 * src1
  echo "fred:"
  echo fred

  echo "calc(fred):"
  echo calc(fred)

  writeEnsemble(fred, "fred.dat")
  echo "\n\nNow try to read"
  let foo = readEnsemble("fred.dat")
  echo "foo:"
  echo foo
  assert(foo =~ fred)
