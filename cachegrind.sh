#!/bin/bash

valgrind --tool=cachegrind ./jacobi 256 1 0.0002
