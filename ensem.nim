## Support for ensemble file format and arithmetic using jackknife/bootstrap propagation of errors.

import
  complex

type
  DataType_t = enum
    RealType = 0,
    ComplexType = 1

  Ensemble_t* = object
    data:      seq[complex]    # Data in either rescaled or raw format
    typ:       DataType_t      # 0 if real and 1 if complex
    nbin:      int             # Number of bins, or data samples, configs, etc
    length:    int             # Number of time-slices of a prop


proc MAX(x: cint; y: cint): cint =
  return if (x < y): y else: x

proc MIN(x: cint; y: cint): cint =
  return if (x < y): x else: y


proc promote_type(src1: DataType_t; src2: DataType_t): DataType_t =
  ## Promote type of data 
  result = RealType
  ## # just to make compiler happy 
  if src1 == RealType and src2 == RealType:
    result = RealType
  elif src1 == RealType and src2 == ComplexType:
    result = ComplexType
  elif src1 == ComplexType and src2 == RealType:
    result = ComplexType
  elif src1 == ComplexType and src2 == ComplexType:
    result = ComplexType
  else:
    quit(stderr, "some unknown types in concatenate")


proc rescale_factor(num: cint): cdouble =
  ## Return the proper rescaling factor 
  ## # Rescaling factor 
  when defined(JACKKNIFE):
    result = - (num - 1)
  else:
    result = sqrt(num - 1)


proc malloc_ensemble(typ: DataType_t; nbin: cint; length: cint): Ensemble_t =
  ## Create a new ensemble of some given number of bins and length 
  if nbin <= 0 or length <= 0:
    quit("invalid input in malloc_ensemble")
  result.data = newSeq[complex](nbin*length)
  if result.data == nil:
    quit("malloc returned NULL in malloc_ensemble")
  result.typ = typ
  result.nbin = nbin
  result.length = length


proc new_ensemble(src: Ensemble_t): Ensemble_t =
  ## Create a new ensemble using parameters from source 
  result = malloc_ensemble(src.typ, src.nbin, src.length)


proc new_len1_ensemble(src: Ensemble_t): Ensemble_t =
  ## # Create a new ensemble of length 1 using the src as input 
  result = malloc_ensemble(src.typ, src.nbin, 1)


proc promote_const_to_ensemble(dval: cdouble; nbin: cint; len: cint): Ensemble_t =
  ## Promote constant to ensemble 
  result = malloc_ensemble(RealType, nbin, len)
  n:cint = 0
  while n < nbin:
    k: cint = 0
    while k < len:
      result.data[k + n * len].real = dval
      result.data[k + n * len].imag = 0.0
      inc(k)
    inc(n)


proc check_two_ensemble(src1: Ensemble_t; src2: Ensemble_t): cint =
  ## # Check if two ensemble have the same parameters 
  if src1.nbin != src2.nbin:
    result = -1
  if src1.length == src2.length:
    result = 0
  elif src1.length == 1 or src2.length > 1:
    result = 1
  elif src1.length > 1 or src2.length == 1:
    result = 2
  else:
    result = -1


proc rescale_ensemble(src: Ensemble_t; factor: cdouble): Ensemble_t =
  ## Return a new rescaled ensemble 
  result = new_ensemble(src)
  var avg: complex_t
  var num: cint = src.nbin
  var len: cint = src.length
  var
    n: cint
    k: cint
  k = 0
  while k < len:
    avg.real = 0.0
    avg.imag = 0.0
    n = 0
    while n < num:
      inc(avg.real, src.data[k + len * n].real)
      inc(avg.imag, src.data[k + len * n].imag)
      inc(n)
    avg.real = avg.real / cast[cdouble](num)
    avg.imag = avg.imag / cast[cdouble](num)
    n = 0
    while n < num:
      result.data[k + len * n].real = avg.real +
          (src.data[k + len * n].real - avg.real) * factor
      result.data[k + len * n].imag = avg.imag +
          (src.data[k + len * n].imag - avg.imag) * factor
      inc(n)
    inc(k)


proc add_const_to_ensemble*(val: cdouble; src2: Ensemble_t): Ensemble_t =
  ## # Add constant on an ensemble 
  result = new_ensemble(src2)
  var num: cint = src2.nbin
  var len: cint = src2.length
  var
    n: cint
    k: cint
  k = 0
  while k < len:
    n = 0
    while n < num:
      result.data[k + len * n].real = val + src2.data[k + len * n].real
      result.data[k + len * n].imag = src2.data[k + len * n].imag
      inc(n)
    inc(k)


proc add_two_ensemble*(src1: Ensemble_t; src2: Ensemble_t): Ensemble_t =
  ## # Add two ensembles 
  var num: cint = src1.nbin
  var len1: cint = src1.length
  var len2: cint = src2.length
  var
    len: cint
    n: cint
    k: cint
    k1: cint
    k2: cint
  case check_two_ensemble(src1, src2)
  of 0, 2:
    result = new_ensemble(src1)
  of 1:
    result = new_ensemble(src2)
  else:
    quit("ensembles not compatible")

  result.typ = promote_type(src1.typ, src2.typ)
  len = MAX(len1, len2)
  k = 0
  while k < len:
    k1 = MIN(k, len1 - 1)
    k2 = MIN(k, len2 - 1)
    n = 0
    while n < num:
      result.data[k + len * n].real = src1.data[k1 + len1 * n].real +
          src2.data[k2 + len2 * n].real
      result.data[k + len * n].imag = src1.data[k1 + len1 * n].imag +
          src2.data[k2 + len2 * n].imag
      inc(n)
    inc(k)


proc subtract_two_ensemble*(src1: Ensemble_t; src2: Ensemble_t): Ensemble_t =
  ## # Subtract two ensembles 
  var num: cint = src1.nbin
  var len1: cint = src1.length
  var len2: cint = src2.length
  var
    len: cint
    n: cint
    k: cint
    k1: cint
    k2: cint
  case check_two_ensemble(src1, src2)
  of 0, 2:
    result = new_ensemble(src1)
  of 1:
    result = new_ensemble(src2)
  else:
    quit("ensembles not compatible")

  result.typ = promote_type(src1.typ, src2.typ)
  len = MAX(len1, len2)
  k = 0
  while k < len:
    k1 = MIN(k, len1 - 1)
    k2 = MIN(k, len2 - 1)
    n = 0
    while n < num:
      result.data[k + len * n].real = src1.data[k1 + len1 * n].real -
          src2.data[k2 + len2 * n].real
      result.data[k + len * n].imag = src1.data[k1 + len1 * n].imag -
          src2.data[k2 + len2 * n].imag
      inc(n)
    inc(k)


proc multiply_const_to_ensemble*(val: cdouble; src2: Ensemble_t): Ensemble_t =
  ## # Multiply a constant on an ensemble 
  result = new_ensemble(src2)
  var num: cint = src2.nbin
  var len: cint = src2.length
  var
    n: cint
    k: cint
  k = 0
  while k < len:
    n = 0
    while n < num:
      result.data[k + len * n].real = val * src2.data[k + len * n].real
      result.data[k + len * n].imag = val * src2.data[k + len * n].imag
      inc(n)
    inc(k)


proc multiply_two_ensemble*(src1: Ensemble_t; src2: Ensemble_t): Ensemble_t =
  ## # Multiply two ensembles 
  var num: cint = src1.nbin
  var len1: cint = src1.length
  var len2: cint = src2.length
  var
    len: cint
    n: cint
    k: cint
    k1: cint
    k2: cint
  case check_two_ensemble(src1, src2)
  of 0, 2:
    result = new_ensemble(src1)
  of 1:
    result = new_ensemble(src2)
  else:
    quit("ensembles not compatible")

  result.typ = promote_type(src1.typ, src2.typ)
  len = MAX(len1, len2)
  k = 0
  while k < len:
    k1 = MIN(k, len1 - 1)
    k2 = MIN(k, len2 - 1)
    n = 0
    while n < num:
      result.data[k + len * n].real = src1.data[k1 + len1 * n].real *
          src2.data[k2 + len2 * n].real -
          src1.data[k1 + len1 * n].imag * src2.data[k2 + len2 * n].imag
      result.data[k + len * n].imag = src1.data[k1 + len1 * n].real *
          src2.data[k2 + len2 * n].imag +
          src1.data[k1 + len1 * n].imag * src2.data[k2 + len2 * n].real
      inc(n)
    inc(k)


proc divide_two_ensemble*(src1: Ensemble_t; src2: Ensemble_t): Ensemble_t =
  ## # Divide two ensembles 
  var num: cint = src1.nbin
  var len1: cint = src1.length
  var len2: cint = src2.length
  var
    len: cint
    n: cint
    k: cint
    k1: cint
    k2: cint
  case check_two_ensemble(src1, src2)
  of 0, 2:
    result = new_ensemble(src1)
  of 1:
    result = new_ensemble(src2)
  else:
    quit("add: ensembles not compatible")

  result.typ = promote_type(src1.typ, src2.typ)
  len = MAX(len1, len2)
  k = 0
  while k < len:
    k1 = MIN(k, len1 - 1)
    k2 = MIN(k, len2 - 1)
    ## # For efficiency, split into real and complex 
    case src2.typ
    of RealType:
      n = 0
      while n < num:
        result.data[k + len * n].real = src1.data[k1 + len1 * n].real div
            src2.data[k2 + len2 * n].real
        result.data[k + len * n].imag = src1.data[k1 + len1 * n].imag div
            src2.data[k2 + len2 * n].real
        inc(n)
    of ComplexType:
      n = 0
      while n < num:
        var denom: cdouble = 1.0 div
            (src2.data[k2 + len2 * n].real * src2.data[k2 + len2 * n].real +
            src2.data[k2 + len2 * n].imag * src2.data[k2 + len2 * n].imag)
        result.data[k + len * n].real = (src1.data[k1 + len1 * n].real *
            src2.data[k2 + len2 * n].real +
            src1.data[k1 + len1 * n].imag * src2.data[k2 + len2 * n].imag) * denom
        result.data[k + len * n].imag = (src1.data[k1 + len1 * n].imag *
            src2.data[k2 + len2 * n].real -
            src1.data[k1 + len1 * n].real * src2.data[k2 + len2 * n].imag) * denom
        inc(n)
    else:
      fprintf(stderr, "something wrong with src2 type: type=%d\x0A", src2.typ)
      exit(1)
    inc(k)


proc negate_ensemble*(src: Ensemble_t): Ensemble_t =
  ## # Negate ensemble 
  result = new_ensemble(src)
  var num: cint = src.nbin
  var len: cint = src.length
  var
    n: cint
    k: cint
  k = 0
  while k < len:
    n = 0
    while n < num:
      result.data[k + len * n].real = - src.data[k + len * n].real
      result.data[k + len * n].imag = - src.data[k + len * n].imag
      inc(n)
    inc(k)


proc real_part_ensemble*(src: Ensemble_t): Ensemble_t =
  ## # Real part ensemble 
  result = new_ensemble(src)
  var num: cint = src.nbin
  var len: cint = src.length
  var
    n: cint
    k: cint
  result.typ = RealType
  k = 0
  while k < len:
    n = 0
    while n < num:
      result.data[k + len * n].real = src.data[k + len * n].real
      result.data[k + len * n].imag = 0.0
      inc(n)
    inc(k)


proc imag_part_ensemble*(src: Ensemble_t): Ensemble_t =
  ## # Imag part ensemble 
  result = new_ensemble(src)
  var num: cint = src.nbin
  var len: cint = src.length
  var
    n: cint
    k: cint
  result.typ = RealType
  k = 0
  while k < len:
    n = 0
    while n < num:
      result.data[k + len * n].real = src.data[k + len * n].imag
      result.data[k + len * n].imag = 0.0
      inc(n)
    inc(k)


proc conj_ensemble*(src: Ensemble_t): Ensemble_t = 
  ## # Conjugate the ensemble 
  result = new_ensemble(src)
  var num: cint = src.nbin
  var len: cint = src.length
  var
    n: cint
    k: cint
  ## # Decide based on the input type 
  case src.typ
  of RealType:
    result.typ = RealType
    k = 0
    while k < len:
      n = 0
      while n < num:
        result.data[k + len * n].real = src.data[k + len * n].real
        result.data[k + len * n].imag = 0.0
        inc(n)
      inc(k)
  of ComplexType:
    result.typ = ComplexType
    k = 0
    while k < len:
      n = 0
      while n < num:
        result.data[k + len * n].real = src.data[k + len * n].real
        result.data[k + len * n].imag = - src.data[k + len * n].imag
        inc(n)
      inc(k)
  else:
    fprintf(stderr, "something wrong with src type: type=%d\x0A", src.typ)
    exit(1)


proc norm2_ensemble*(src: Ensemble_t): Ensemble_t =
  ## # Norm2 the ensemble 
  result = new_ensemble(src)
  var num: cint = src.nbin
  var len: cint = src.length
  var
    n: cint
    k: cint
  ## # Decide based on the input type 
  case src.typ
  of RealType:
    result.typ = RealType
    k = 0
    while k < len:
      n = 0
      while n < num:
        var re: cdouble = src.data[k + len * n].real
        result.data[k + len * n].real = re * re
        result.data[k + len * n].imag = 0.0
        inc(n)
      inc(k)
  of ComplexType:
    result.typ = RealType
    k = 0
    while k < len:
      n = 0
      while n < num:
        var re: cdouble = src.data[k + len * n].real
        var im: cdouble = src.data[k + len * n].imag
        result.data[k + len * n].real = re * re + im * im
        result.data[k + len * n].imag = 0.0
        inc(n)
      inc(k)
  else:
    quit("something wrong with src type: type= " & $src.typ)


proc atan2_ensemble*(src1: Ensemble_t; src2: Ensemble_t): Ensemble_t =
  ## # Call atan2 on two real ensembles 
  var num: cint = src1.nbin
  var len1: cint = src1.length
  var len2: cint = src2.length
  var
    len: cint
    n: cint
    k: cint
    k1: cint
    k2: cint
  case check_two_ensemble(src1, src2)
  of 0, 2:
    result = new_ensemble(src1)
  of 1:
    result = new_ensemble(src2)
  else:
    quit("cmplx: ensembles not compatible")

  if not (src1.typ == RealType and src2.typ == RealType):
    quit("atan2 requires both ensembles real")

  result.typ = RealType
  len = MAX(len1, len2)
  k = 0
  while k < len:
    k1 = MIN(k, len1 - 1)
    k2 = MIN(k, len2 - 1)
    n = 0
    while n < num:
      result.data[k + len * n].real = atan2(src1.data[k1 + len1 * n].real,
                                   src2.data[k2 + len2 * n].real)
      result.data[k + len * n].imag = 0.0
      inc(n)
    inc(k)qu


proc cmplx_ensemble*(src1: Ensemble_t; src2: Ensemble_t): Ensemble_t =
  ## # Build complex from two real ensembles 
  var num: cint = src1.nbin
  var len1: cint = src1.length
  var len2: cint = src2.length
  var
    len: cint
    n: cint
    k: cint
    k1: cint
    k2: cint
  case check_two_ensemble(src1, src2)
  of 0, 2:
    result = new_ensemble(src1)
  of 1:
    result = new_ensemble(src2)
  else:
    quit("cmplx: ensembles not compatible")

  if not (src1.typ == RealType and src2.typ == RealType):
    quit("cmplx requires both ensembles real")

  result.typ = ComplexType
  len = MAX(len1, len2)
  k = 0
  while k < len:
    k1 = MIN(k, len1 - 1)
    k2 = MIN(k, len2 - 1)
    n = 0
    while n < num:
      result.data[k + len * n].real = src1.data[k1 + len1 * n].real
      result.data[k + len * n].imag = src2.data[k2 + len2 * n].real
      inc(n)
    inc(k)


proc timesI_ensemble*(src: Ensemble_t): Ensemble_t =
  ## # Multiply ensemble by I 
  result = new_ensemble(src)
  var num: cint = src.nbin
  var len: cint = src.length
  var
    n: cint
    k: cint
  result.typ = ComplexType
  case src.typ
  of RealType:
    k = 0
    while k < len:
      n = 0
      while n < num:
        result.data[k + len * n].real = 0.0
        result.data[k + len * n].imag = src.data[k + len * n].real
        inc(n)
      inc(k)
  of ComplexType:
    k = 0
    while k < len:
      n = 0
      while n < num:
        result.data[k + len * n].real = - src.data[k + len * n].imag
        result.data[k + len * n].imag = src.data[k + len * n].real
        inc(n)
      inc(k)
  else:
    quit("something wrong with ensemble: type= " & $src.typ)


proc read_fileptr_ensemble*(fp: ptr FILE; name: string): Ensemble_t =
  ## # Read ensemble 
  var num: cint
  var len: cint
  var
    n: cint
    k: cint
  var
    junk: cint
    ncol: cint
    ttype: cint
  var typ: DataType_t
  const
    MAXLINE = 1000
  var crap: array[MAXLINE, char]
  n = 0
  while (crap[inc(n)] = fgetc(fp)) != '\x0A':
    if n == MAXLINE:
      fprintf(stderr, "Header line too long\x0A")
      exit(1)
  if sscanf(crap, "%d %d %d %d %d", addr(num), addr(len), addr(ttype), addr(junk),
           addr(ncol)) != 5:
    fprintf(stderr, "error reading parameters from %s\x0A", name)
    exit(1)
  typ = ttype
  if ncol != 1:
    fprintf(stderr, "only support 1 column in %s\x0A", name)
    exit(1)
  if typ != RealType and typ != ComplexType:
    fprintf(stderr, "error in type for file %s\x0A", name)
    exit(1)
  result = malloc_ensemble(typ, num, len)
  n = 0
  while n < num:
    var t: cint
    case typ
    of RealType:
      k = 0
      while k < len:
        if fscanf(fp, "%d %lf", addr(t), addr((result.data[k + len * n].real))) != 2:
          fprintf(stderr, "error reading data from %s\x0A", name)
          exit(1)
        if k != t:
          fprintf(stderr, "error reading time slice data from %s\x0A", name)
          exit(1)
        result.data[k + len * n].imag = 0.0
        inc(k)
    of ComplexType:
      k = 0
      while k < len:
        if fscanf(fp, "%d %lf %lf", addr(t), addr((result.data[k + len * n].real)),
                 addr((result.data[k + len * n].imag))) != 3:
          fprintf(stderr, "error reading data from %s\x0A", name)
          exit(1)
        if k != t:
          fprintf(stderr, "error reading time slice data from %s\x0A", name)
          exit(1)
        inc(k)
    else:
      fprintf(stderr, "something wrong with ensemble: type=%d\x0A", typ)
      exit(1)
    inc(n)


proc read_ensemble*(name: string): Ensemble_t =
  ## # Read ensemble from named file 
  var fp: FILE
  fp = fopen(name, "r")
  if fp == nil:
    fprintf(stderr, "file %s does not exist\x0A", name)
    exit(1)
  result = read_fileptr_ensemble(fp, name)
  fclose(fp)


proc write_ensemble*(name: string; src: Ensemble_t) =
  ## # Write ensemble 
  var typ: DataType_t = src.typ
  var num: cint = src.nbin
  var len: cint = src.length
  var
    n: cint
    k: cint
  var fp: ptr FILE
  fp = fopen(name, "w")
  if fp == nil:
    fprintf(stderr, "error opening file %s\x0A", name)
    exit(1)
  fprintf(fp, "%d %d %d 0 1\x0A", num, len, typ)
  n = 0
  while n < num:
    case typ
    of RealType:
      k = 0
      while k < len:
        fprintf(fp, "%d %.12g\x0A", k, src.data[k + len * n].real)
        inc(k)
    of ComplexType:
      k = 0
      while k < len:
        fprintf(fp, "%d %.12g %.12g\x0A", k, src.data[k + len * n].real,
                src.data[k + len * n].imag)
        inc(k)
    else:
      fprintf(stderr, "something wrong with ensemble: type=%d\x0A", typ)
      exit(1)
    inc(n)
  fclose(fp)


proc apply_func_ensemble(funcptr: proc (): cdouble; src: Ensemble_t): Ensemble_t =
  ## # Apply function to ensemble 
  result = new_ensemble(src)
  var num: cint = src.nbin
  var len: cint = src.length
  var
    n: cint
    k: cint
  if src.typ != RealType:
    quit("only func(real) not supported")
  k = 0
  while k < len:
    n = 0
    while n < num:
      result.data[k + len * n].real = (funcptr)(src.data[k + len * n].real)
      result.data[k + len * n].imag = 0.0
      inc(n)
    inc(k)


proc apply_pow_ensemble(src: Ensemble_t; val: cdouble): Ensemble_t =
  ## # Apply pow(src,const) to ensemble 
  result = new_ensemble(src)
  var num: cint = src.nbin
  var len: cint = src.length
  var
    n: cint
    k: cint
  if src.typ != RealType:
    quit("only func(real) not supported")
  k = 0
  while k < len:
    n = 0
    while n < num:
      result.data[k + len * n].real = pow(src.data[k + len * n].real, val)
      result.data[k + len * n].imag = 0.0
      inc(n)
    inc(k)


proc calc_real_ensemble(src: Ensemble_t) =
  ## # Calculate mean, err and err/mean for data 
  var avg: cdouble
  var err: cdouble
  var rat: cdouble
  var diff: cdouble
  var num: cint = src.nbin
  var len: cint = src.length
  var
    n: cint
    k: cint
  var fmt: ptr cstring = ["%d   %g %g   %g\x0A",
                     "%4d   %- 13.6g %- 13.6g   %- 13.6g\x0A"]
  var dofmt: cint
  dofmt = if (len > 1): 1 else: 0
  k = 0
  while k < len:
    avg = 0.0
    err = 0.0
    n = 0
    while n < num:
      inc(avg, src.data[k + len * n].real)
      inc(n)
    avg = avg / cast[cdouble](num)
    n = 0
    while n < num:
      diff = (src.data[k + len * n].real - avg)
      inc(err, diff * diff)
      inc(n)
    err = sqrt(err div (double)((num - 1) * num))
    if avg != 0.0: rat = err div avg
    else: rat = 0.0
    printf(fmt[dofmt], k, avg, err, rat)
    inc(k)


proc calc_complex_ensemble*(src: Ensemble_t) =
  ## # Calculate mean, err and err/mean for data 
  var avg: complex_t
  var err: cdouble
  var diff: complex_t
  var num: cint = src.nbin
  var len: cint = src.length
  var
    n: cint
    k: cint
  var fmt: ptr cstring = ["%d   ( %g , %g )   %g\x0A",
                     "%4d   ( %- 13.6g , %- 13.6g )   %- 13.6g\x0A"]
  var dofmt: cint
  dofmt = if (len > 1): 1 else: 0
  k = 0
  while k < len:
    avg.real = avg.imag = 0.0
    err = 0.0
    n = 0
    while n < num:
      inc(avg.real, src.data[k + len * n].real)
      inc(avg.imag, src.data[k + len * n].imag)
      inc(n)
    avg.real = avg.real / cast[cdouble](num)
    avg.imag = avg.imag / cast[cdouble](num)
    n = 0
    while n < num:
      diff.real = (src.data[k + len * n].real - avg.real)
      diff.imag = (src.data[k + len * n].imag - avg.imag)
      inc(err, diff.real * diff.real + diff.imag * diff.imag)
      inc(n)
    err = sqrt(err div (double)((num - 1) * num))
    printf(fmt[dofmt], k, avg.real, avg.imag, err)
    inc(k)


proc calc_ensemble*(src: Ensemble_t) =
  ## # Calculate mean, err and err/mean for data 
  case src.typ
  of RealType:
    calc_real_ensemble(src)
  of ComplexType:
    calc_complex_ensemble(src)
  else:
    fprintf(stderr, "something wrong with ensemble: type=%d\x0A", src.typ)
    exit(1)


proc print_ensemble*(src: Ensemble_t) =
  ## # Print the ensemble (maybe for a pipe??) 
  var typ: DataType_t = src.typ
  var num: cint = src.nbin
  var len: cint = src.length
  var
    n: cint
    k: cint
  var fp: ptr FILE
  fp = stdout
  if fp == nil:
    fprintf(stderr, "something wrong with stdout\x0A")
    exit(1)
  fprintf(fp, "%d %d %d 0 1\x0A", num, len, typ)
  case typ
  of RealType:
    n = 0
    while n < num:
      k = 0
      while k < len:
        fprintf(fp, "%d %.12g\x0A", k, src.data[k + len * n].real)
        inc(k)
      inc(n)
  of ComplexType:
    n = 0
    while n < num:
      k = 0
      while k < len:
        fprintf(fp, "%d %.12g %.12g\x0A", k, src.data[k + len * n].real,
                src.data[k + len * n].imag)
        inc(k)
      inc(n)
  else:
    fprintf(stderr, "something wrong with ensemble: type=%d\x0A", typ)
    exit(1)


proc shift_ensemble*(src: Ensemble_t; sh: cint): Ensemble_t =
  ## # Shift an ensemble in some direction dropping bits from the end 
  result = new_ensemble(src)
  var num: cint = src.nbin
  var len: cint = src.length
  var
    n: cint
    k: cint
    kk: cint
  if abs(sh) > len:
    fprintf(stderr, "shift: do not allow the shift greater than the length\x0A")
    exit(1)
  n = 0
  while n < num:
    k = 0
    while k < len:
      kk = (k + len + sh) mod len
      result.data[k + len * n].real = src.data[kk + len * n].real
      result.data[k + len * n].imag = src.data[kk + len * n].imag
      inc(k)
    ## # Clean out the ends 
    if sh > 0:
      k = len - sh
      while k < len:
        result.data[k + len * n].real = 0.0
        result.data[k + len * n].imag = 0.0
        inc(k)
    elif sh < 0:
      k = 0
      while k < - sh:
        result.data[k + len * n].real = 0.0
        result.data[k + len * n].imag = 0.0
        inc(k)
    inc(n)


proc cshift_ensemble*(src: Ensemble_t; sh: cint): Ensemble_t =
  ## # Periodic (circular) shift an ensemble in some direction 
  result = new_ensemble(src)
  var num: cint = src.nbin
  var len: cint = src.length
  var
    n: cint
    k: cint
    kk: cint
  n = 0
  while n < num:
    k = 0
    while k < len:
      kk = (k + len + sh) mod len
      result.data[k + len * n].real = src.data[kk + len * n].real
      result.data[k + len * n].imag = src.data[kk + len * n].imag
      inc(k)
    inc(n)


proc extract_ensemble*(src: Ensemble_t; elem_i: cint; elem_f: cint): Ensemble_t =
  ## # Extract a range of time slices from an ensemble 
  result.typ = src.typ
  var num: cint = src.nbin
  var len: cint = src.length
  var dlen: cint = elem_f - elem_i
  var
    k: cint
    n: cint
  if elem_i < 0 or elem_i >= len:
    fprintf(stderr, "index element out of bounds of ensemble: %d\x0A", elem_i)
    exit(1)
  if elem_f < 0 or elem_f >= len:
    fprintf(stderr, "index element out of bounds of ensemble: %d\x0A", elem_f)
    exit(1)
  if elem_f < elem_i:
    fprintf(stderr, "index element out of order in ensemble: %d %d\x0A", elem_i,
            elem_f)
    exit(1)
  dlen = elem_f - elem_i + 1
  result = malloc_ensemble(typ, num, dlen)
  n = 0
  while n < num:
    k = 0
    while k < dlen:
      result.data[k + dlen * n].real = src.data[k + elem_i + len * n].real
      result.data[k + dlen * n].imag = src.data[k + elem_i + len * n].imag
      inc(k)
    inc(n)


proc concatenate_ensemble*(src1: Ensemble_t; src2: Ensemble_t): Ensemble_t =
  ## # Concatenate two ensembles 
  var
    k: cint
    n: cint
  var
    len: cint
    num: cint
  var typ: DataType_t
  if src1.nbin != src2.nbin:
    quit("Ensembles not compatible for concatenation")

  len = src1.length + src2.length
  num = src1.nbin
  typ = promote_type(src1.typ, src2.typ)
  result = malloc_ensemble(typ, num, len)
  n = 0
  while n < num:
    k = 0
    while k < src1.length:
      result.data[k + n * len].real = src1.data[k + n * src1.length].real
      result.data[k + n * len].imag = src1.data[k + n * src1.length].imag
      inc(k)
    k = 0
    while k < src2.length:
      result.data[k + src1.length + n * len].real = src2.data[k + n * src2.length].real
      result.data[k + src1.length + n * len].imag = src2.data[k + n * src2.length].imag
      inc(k)
    inc(n)
