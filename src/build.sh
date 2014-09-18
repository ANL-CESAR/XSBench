#!/bin/sh
java -classpath /home/alund/openarc/lib/cetus.jar:/home/alund/openarc/lib/antlr.jar openacc.exec.ACC2GPUDriver -macro=OPENACC -verbosity=2 *.c
