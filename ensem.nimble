# Package
version       = "1.3.0"
author        = "Robert Edwards"
description   = "Support for ensemble file format and arithmetic using jackknife/bootstrap propagation of errors"
license       = "BSD3"
srcDir        = "src"

# Dependencies
requires "nim >= 0.19.0"

# Builds
task test, "Run the test suite":
  exec "nim c -r tests/test_ensem"

task docgen, "Regenerate the documentation":
  exec "nim doc2 --out:docs/ensem.html src/ensem.nim"


