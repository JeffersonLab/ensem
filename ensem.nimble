# Package

version       = "1.2.0"
author        = "Robert Edwards"
description   = "Support for ensemble file format and arithmetic using jackknife/bootstrap propagation of errors"
license       = "MIT"

# Dependencies

requires "nim >= 0.17.0"

# Builds
skipDirs = @["tests"]

task test, "Runs the test suite":
  exec "nim c -r tests/test_ensem"


