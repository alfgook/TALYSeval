# TALYSeval - R package

This package provides functionality to set up [TALYS](www.talys.eu) calculations
and retrieve their results.

## Installation

Run the following commands in a terminal:
```
git clone https://github.com/gschnabel/TALYSeval.git
R CMD INSTALL TALYSeval
```

## Basic usage

In order to set up a calculation, one defines a named list
in which the keys are the TALYS keywords and the 
associated values the parameter values:

```
library(TALYSeval)

parameterSet <- list(
    projectile = "n",
    element = "Fe",
    mass = 56,
    energy = c(1,3,5),
    endf = "y"
)
```
Note that the energy keyword takes a numeric vector
with the energies instead of a filename.

With the defined list of TALYS parameters, 
the calculation can be set up in the directory, e.g.,
`testcalcdir` which must exist before executing:

```
talysModel <- createModelTALYS()
talysModel$prepare("testcalcdir", parameterSet)
```
TALYS can be run from R using the command:
```
setwd("testcalcdir")
system("talys < input > output")
```

After the successful calculation, predictions can be retrieved.
First, a datatable with the desired result types needs to be defined:
```
library(data.table)
resultSpec <- data.table(REAC="CS/EL", L1=seq(1,3,length=10), L2=0, L3=0)
```

The results of the calculation can be retrieved by invoking:
```
result <- talysModel$read("testcalc", resultSpec, packed=FALSE)
print(result)
```

