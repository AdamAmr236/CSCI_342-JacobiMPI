#!/bin/bash

objdump \
  --disassemble \
  --demangle \
  --no-show-raw-insn \
  obj/worker.o > obj/worker.s
