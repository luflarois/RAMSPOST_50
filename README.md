# RAMSPOST_5.0

System to pos processing BRAMS' model analisys output

### Authors

by Ricardo Almeida de Siqueira, Madeleine Sánchez Gácita - Jan 2012

Improvements by Saulo Ribeiro de Freitas - Set 2015

Improvements by Luiz Flavio Rodrigues - Ago 2025


### Abstract

O RAMSPOST (RAMS-POST processing) é um pacote para a reformatação da saída do
modelo BRAMS visando à geração de gráficos de variáveis ambientais. O BRAMS é um
modelo regional derivado do RAMS (http://brams.cptec.inpe.br/) e seus arquivos de saída
possuem diferentes informações atmosféricas (como os componentes do vetor velocidade do
vento, pressão atmosférica e temperatura), informações geográficas (como a topografia) e
outros.
Uma vez que o BRAMS termina uma rodada, um conjunto de arquivos chamado de "análises"
é gerado num formato de arquivo chamado vfm. O vfm é um formato de arquivo que serve para
a geração dos arquivos de entrada e saída do BRAMS onde, no caso das análises, estão as
informações de todas as variáveis de saída retornadas pelo modelo. O RAMSPOST usa os
arquivos de análises vfm como entrada para produzir os arquivos de saída correspondentes ao
formato de entrada do software GrADS (Grid Analysis and Display System) para a posterior
visualização gráfica. O GrADS é uma ferramenta que permite o acesso, a manipulação e a
visualização de dados de ciências da Terra (http://grads.iges.org/grads/grads.html).

### Estrutura

O RAMSPOST é composto pelos seguintes arquivos:
* Makefile_50
* src/anheader.f90
* include_ramspost.mk
* ramspost.inp
* src/ramspost_A.f90
* src/ramspost_B.f90
* src/ramspost_C.f90
* src/ramspost_D.f90
e pelos diretórios:
* LIB
* include
* util

## How to compile/build

This file contains instructions for compiling RAMSPOST version 5.0.
The goal is pos-processing BRAMS 5.1 output files (analysis) for visualization with GrADS software.

Follow the steps:
1. making  library "libutils-ramspost.a"
  a) go to "LIB" directory
  b) set your fortran and C compilers using the file named "include.mk.opt.ramspost"
  c) execute the command:  make -f Make.utils.opt
2. making the executable "ramspost_50"
  a) set your fortran and C compilers using the file named "include_ramspost.mk"
  b) execute the command:  make -f Makefile_50
  
The file "ramspost.inp" is a namelist read by the executable "ramspost_50" which will create 
GrADS files with the variables defined in the namelist.
For a complete list of the variables avalaible, look at the file 
"README_var_library_for_BRAMS5.1"



