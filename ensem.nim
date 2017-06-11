## Support for ensemble file format and arithmetic using jackknife/bootstrap propagation of errors.

import
  complex, math, streams, parseutils, strutils

type
  DataType_t* = enum
    RealType = 0,
    ComplexType = 1

  EnsemType_t* = enum
    EnsemJackknife = 10,
    EnsemBootstrap = 11

  Ensemble_t* = object
    data:        seq[Complex]    # Data in either rescaled or raw format
    typ:         DataType_t      # 0 if real and 1 if complex
    nbin:        int             # Number of bins, or data samples, configs, etc
    Lt:          int             # Number of time-slices of a prop
#    ens:         EnsemType_t     # Type of fluctuations for the ensemble


proc norm2(x: Complex): float =
  result = x.re * x.re + x.im * x.im

proc MAX(x: int; y: int): int =
  return if (x < y): y else: x

proc MIN(x: int; y: int): int =
  return if (x < y): x else: y


proc numElem*(src: Ensemble_t): int = 
  ## Get the time-extent
  result = src.Lt


proc size*(src: Ensemble_t): int = 
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


proc newRealEnsemble*(val: float; nbin: int; len: int): Ensemble_t =
  ## Promote constant to ensemble 
  result = newEnsemble(RealType, nbin, len)
  var n = 0
  while n < nbin:
    var k = 0
    while k < len:
      result.data[k + n * len].re = val
      result.data[k + n * len].im = 0.0
      inc(k)
    inc(n)


proc newComplexEnsemble*(val: Complex; nbin: int; len: int): Ensemble_t =
  ## Promote constant to ensemble 
  result = newEnsemble(ComplexType, nbin, len)
  var n = 0
  while n < nbin:
    var k = 0
    while k < len:
      result.data[k + n * len] = val
      inc(k)
    inc(n)


proc check_two_ensemble(src1: Ensemble_t; src2: Ensemble_t): int =
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


proc rescale_ensemble(src: Ensemble_t; factor: float): Ensemble_t =
  ## Return a new rescaled ensemble 
  result = newEnsemble(src)
  var avg: Complex
  var num: int = src.nbin
  var Lt:  int = src.Lt

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
  


proc `+`*(val: float; src2: Ensemble_t): Ensemble_t =
  ## Add constant on an ensemble 
  result = newEnsemble(src2)
  var Lt:  int = src2.Lt
  for k in 0..Lt-1:
    for n in 0..src2.nbin-1:
      result.data[k + Lt * n] = val + src2.data[k + Lt * n]


proc `+`*(val: Complex; src2: Ensemble_t): Ensemble_t =
  ## Add constant on an ensemble 
  result = newEnsemble(src2)
  var Lt:  int = src2.Lt
  for k in 0..Lt-1:
    for n in 0..src2.nbin-1:
      result.data[k + Lt * n] = val + src2.data[k + Lt * n]


proc `+`*(src1: Ensemble_t; val: float): Ensemble_t =
  ## Add constant on an ensemble 
  result = val + src1


proc `+`*(src1: Ensemble_t; val: Complex): Ensemble_t =
  ## Add constant on an ensemble 
  result = val + src1


proc `+`*(src1: Ensemble_t; src2: Ensemble_t): Ensemble_t =
  ## Add two ensembles 
  var num: int = src1.nbin
  var Lt1: int = src1.Lt
  var Lt2: int = src2.Lt

  case check_two_ensemble(src1, src2)
  of 0, 2:
    result = newEnsemble(src1)
  of 1:
    result = newEnsemble(src2)
  else:
    quit("ensembles not compatible")

  result.typ = promote_type(src1.typ, src2.typ)
  var Lt = MAX(Lt1, Lt2)
  for k in 0..Lt-1:
    var k1 = MIN(k, Lt1 - 1)
    var k2 = MIN(k, Lt2 - 1)
    for n in 0..num-1:
      result.data[k + Lt * n] = src1.data[k1 + Lt1 * n] + src2.data[k2 + Lt2 * n]


proc `-`*(src1: Ensemble_t; src2: Ensemble_t): Ensemble_t =
  ## Subtract two ensembles 
  var num: int = src1.nbin
  var Lt1: int = src1.Lt
  var Lt2: int = src2.Lt
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


proc `*`*(val: float; src2: Ensemble_t): Ensemble_t =
  ## Multiply a constant on an ensemble 
  result = newEnsemble(src2)
  var Lt:  int = src2.Lt
  for k in 0..Lt-1:
    for n in 0..src2.nbin-1:
      result.data[k + Lt * n] = val * src2.data[k + Lt * n]


proc `*`*(val: Complex; src2: Ensemble_t): Ensemble_t =
  ## Multiply a constant on an ensemble 
  result = newEnsemble(src2)
  var Lt:  int = src2.Lt
  for k in 0..Lt-1:
    for n in 0..src2.nbin-1:
      result.data[k + Lt * n] = val * src2.data[k + Lt * n]


proc `*`*(src1: Ensemble_t; val: float): Ensemble_t =
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
  var num: int = src1.nbin
  var Lt1: int = src1.Lt
  var Lt2: int = src2.Lt
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
  var num: int = src1.nbin
  var Lt1: int = src1.Lt
  var Lt2: int = src2.Lt

  case check_two_ensemble(src1, src2)
  of 0, 2:
    dest = newEnsemble(src1)
  of 1:
    dest = newEnsemble(src2)
  else:
    quit("add: ensembles not compatible")

  dest.typ = promote_type(src1.typ, src2.typ)
  var Lt = MAX(Lt1, Lt2)

  for k in 0..Lt-1:
    var k1 = MIN(k, Lt1 - 1)
    var k2 = MIN(k, Lt2 - 1)
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
        var denom: float = 1.0 /
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
  var Lt:  int = src.Lt
  for k in 0..Lt-1:
    for n in 0..src.nbin-1:
      result.data[k + Lt * n] = - src.data[k + Lt * n]


proc real*(src: Ensemble_t): Ensemble_t =
  ## Real part ensemble 
  result = newEnsemble(src)
  var Lt:  int = src.Lt
  result.typ = RealType
  for k in 0..Lt-1:
    for n in 0..src.nbin-1:
      result.data[k + Lt * n].re = src.data[k + Lt * n].re
      result.data[k + Lt * n].im = 0.0


proc imag*(src: Ensemble_t): Ensemble_t =
  ## Imag part ensemble 
  result = newEnsemble(src)
  var Lt:  int = src.Lt
  result.typ = RealType
  for k in 0..Lt-1:
    for n in 0..src.nbin-1:
      result.data[k + Lt * n].re = src.data[k + Lt * n].im
      result.data[k + Lt * n].im = 0.0


proc conjugate*(src: Ensemble_t): Ensemble_t = 
  ## Conjugate the ensemble 
  result = newEnsemble(src)
  var Lt:  int = src.Lt

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
  var Lt:  int = src.Lt

  ## Decide based on the input type 
  case src.typ
  of RealType:
    dest.typ = RealType
    for k in 0..Lt-1:
      for n in 0..src.nbin-1:
        var re: float = src.data[k + Lt * n].re
        dest.data[k + Lt * n].re = re * re
        dest.data[k + Lt * n].im = 0.0

  of ComplexType:
    dest.typ = RealType
    for k in 0..Lt-1:
      for n in 0..src.nbin-1:
        var re: float = src.data[k + Lt * n].re
        var im: float = src.data[k + Lt * n].im
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
  var num: int = src1.nbin
  var Lt1: int = src1.Lt
  var Lt2: int = src2.Lt

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
  var Lt = MAX(Lt1, Lt2)
  for k in 0..Lt-1:
    var k1 = MIN(k, Lt1 - 1)
    var k2 = MIN(k, Lt2 - 1)
    for n in 0..num-1:
      result.data[k + Lt * n].re = arctan2(src1.data[k1 + Lt1 * n].re,
                                   src2.data[k2 + Lt2 * n].re)
      result.data[k + Lt * n].im = 0.0


proc cmplx*(src1: Ensemble_t; src2: Ensemble_t): Ensemble_t =
  ## Build complex from two real ensembles 
  var num: int = src1.nbin
  var Lt1: int = src1.Lt
  var Lt2: int = src2.Lt

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
  var Lt = MAX(Lt1, Lt2)
  for k in 0..Lt-1:
    var k1 = MIN(k, Lt1 - 1)
    var k2 = MIN(k, Lt2 - 1)
    for n in 0..num-1:
      result.data[k + Lt * n].re = src1.data[k1 + Lt1 * n].re
      result.data[k + Lt * n].im = src2.data[k2 + Lt2 * n].re


proc timesI*(src: Ensemble_t): Ensemble_t =
  ## Multiply ensemble by I 
  result = newEnsemble(src)
  var num: int = src.nbin
  var Lt:  int = src.Lt

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
  var fs = newFileStream(name, fmRead)
  
  if isNil(fs):
    quit("Error opening file = " & name)

  var line: string
  if not fs.readLine(line):
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


proc writeEnsemble*(name: string; src: Ensemble_t) =
  ## Write ensemble 
  var typ: DataType_t = src.typ
  var num = src.nbin
  var Lt  = src.Lt
  var
    n: int
    k: int
  var fp: FILE
  if not open(fp, name):
    quit("error opening file " & name)

  writeLine(fp, num, Lt, typ, "0 1")
  n = 0
  while n < num:
    case typ
    of RealType:
      k = 0
      while k < Lt:
        writeLine(fp, k, src.data[k + Lt * n].re)
        inc(k)
    of ComplexType:
      k = 0
      while k < Lt:
        writeLine(fp, k, src.data[k + Lt * n].re, src.data[k + Lt * n].im)
        inc(k)
    else:
      quit("something wrong with ensemble: type= " & $typ)
    inc(n)
  close(fp)


#
#proc apply_func_ensemble(funcptr: proc (): float; src: Ensemble_t): Ensemble_t =
#  ## Apply function to ensemble 
#  result = newEnsemble(src)
#  var num = src.nbin
#  var Lt  = src.Lt
#  var
#    n: int
#    k: int
#  if src.typ != RealType:
#    quit("only func(real) not supported")
#  k = 0
#  while k < Lt:
#    n = 0
#    while n < num:
#      result.data[k + Lt * n].re = (funcptr)(src.data[k + Lt * n].re)
#      result.data[k + Lt * n].im = 0.0
#      inc(n)
#    inc(k)
#

proc apply_pow_ensemble(src: Ensemble_t; val: float): Ensemble_t =
  ## Apply pow(src,const) to ensemble 
  result = newEnsemble(src)
  var num: int = src.nbin
  var Lt:  int = src.Lt
  var
    n: int
    k: int
  if src.typ != RealType:
    quit("only func(real) not supported")
  k = 0
  while k < Lt:
    n = 0
    while n < num:
      result.data[k + Lt * n].re = pow(src.data[k + Lt * n].re, val)
      result.data[k + Lt * n].im = 0.0
      inc(n)
    inc(k)


proc calc_real_ensemble(src: Ensemble_t) =
  ## Calculate mean, err and err/mean for data 
  var avg: float
  var err: float
  var rat: float
  var diff: float
  var num: int = src.nbin
  var Lt:  int = src.Lt
#  let fmt = ["%d   %g %g   %g", "%4d   %- 13.6g %- 13.6g   %- 13.6g"]
#  var dofmt: int
#  dofmt = if (Lt > 1): 1 else: 0
  for k in 0..Lt-1:
    avg = 0.0
    err = 0.0
    for n in 0..num-1:
      avg += src.data[k + Lt * n].re

    avg /= float(num)
    for n in 0..num-1:
      diff = (src.data[k + Lt * n].re - avg)
      err += diff * diff

    err = sqrt(err / float((num - 1) * num))
    if avg != 0.0: rat = err / avg
    else: rat = 0.0
    echo k, "   ", avg, " ", err, "   ", rat


proc calc_complex_ensemble(src: Ensemble_t) =
  ## Calculate mean, err and err/mean for data 
  var avg: Complex
  var err: float
  var diff: Complex
  var num: int = src.nbin
  var Lt:  int = src.Lt
#  var fmt: ptr cstring = ["%d   ( %g , %g )   %g",
#                     "%4d   ( %- 13.6g , %- 13.6g )   %- 13.6g"]
#  var dofmt: int
#  dofmt = if (Lt > 1): 1 else: 0
  for k in 0..Lt-1:
    avg = (0.0, 0.0)
    err = 0.0
    for n in 0..num-1:
      avg += src.data[k + Lt * n]

    avg /= float(num)
    for n in 0..num-1:
      diff = src.data[k + Lt * n] - avg
      err += norm2(diff)

    err = sqrt(err / float((num - 1) * num));
    echo k, "   (", avg.re, " , ", avg.im, ")   ", err;


proc calc*(src: Ensemble_t) =
  ## Calculate mean, err and err/mean for data 
  case src.typ
  of RealType:
    calc_real_ensemble(src)
  of ComplexType:
    calc_complex_ensemble(src)
  else:
    quit("something wrong with ensemble: type= " & $src.typ)


proc print*(src: Ensemble_t) =
  ## Print the ensemble (maybe for a pipe??) 
  var typ: DataType_t = src.typ
  var num: int = src.nbin
  var Lt:  int = src.Lt

#  writeLine(fp, "%d %d %d 0 1", num, Lt, typ)
#  writeLine(stdout, "%d %d %d 0 1", num, Lt, typ)
  echo num, " ", Lt, " ", int(typ), " ", 0, " ", 1
  case typ
  of RealType:
    for n in 0..num-1:
      for k in 0..Lt-1:
#        writeLine(stdout, "%d %.12g", k, src.data[k + Lt * n].re)
        echo k, " ", src.data[k + Lt * n].re

  of ComplexType:
    for n in 0..num-1:
      for k in 0..Lt-1:
#        writeLine(stdout, "%d %.12g %.12g", k, src.data[k + Lt * n].re, src.data[k + Lt * n].im)
        echo k, " ", src.data[k + Lt * n].re, " ", src.data[k + Lt * n].im

  else:
    quit("something wrong with ensemble: type= " & $typ)


proc shift*(src: Ensemble_t; sh: int): Ensemble_t =
  ## Shift an ensemble in some direction dropping bits from the end 
  result = newEnsemble(src)
  var num: int = src.nbin
  var Lt:  int = src.Lt
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
  result = newEnsemble(src)
  var num: int = src.nbin
  var Lt:  int = src.Lt
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


proc extract*(src: Ensemble_t; elem_i: int; elem_f: int): Ensemble_t =
  ## Extract a range of time slices from an ensemble 
  result.typ = src.typ
  result.typ = src.typ
  var num: int = src.nbin
  var Lt:  int = src.Lt
  var dLt: int = elem_f - elem_i
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
  let src1 = newRealEnsemble(5.0, 3, 4)
  let src2 = newRealEnsemble(17.0, 3, 4)
  let src3 = newComplexEnsemble((5.0,7.3), 3, 4)

#  echo "src1:"
#  print(src1)
  echo "src2:"
  print(src2)
  echo "src3:"
  print(src3)

  var fred = src3 * src2
  echo "fred:"
  print(fred)

  echo "calc(fred):"
  calc(fred)
