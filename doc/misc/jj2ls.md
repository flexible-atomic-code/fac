# Using jj2LSJ with FAC

`BasisTable` function can produce GRASP2K files in order to make use of its `jj2LSJ` module.

Here is an example usage of this function,

```
SetAtom('Be')

Config('g1', '1*2 2*2')
Config('g3', '1*2 2*1 3*1')
Config('g4', '1*2 2*1 4*1')

ConfigEnergy(0)
OptimizeRadial('g1')
ConfigEnergy(1)

Structure('Be.en', ['g1']
Structure('Be.en', ['g3']
Structure('Be.en', ['g4']

BasisTable('Be', 10)
```

which generates `Be_001.c, Be_002.c, ..., Be_008.c` and `Be_001.cm, Be_002.cm, ..., Be_008.cm`.

Here, the suffix `_00*` indicates number of blocks in the output. Since GRASP2K can handle only one block in each run, we need to pass these files separately.

For example, the first few lines of `Be_002.c` is something like
```
Core subshells:

Peel subshells:
  1s   2s   2p-  2p   4s   4p-  4p   4d-  4d   4f-  4f
CSF(s):
  1s ( 2)  2p ( 1)  4p ( 1)
               3/2      3/2
                           0+
  1s ( 2)  2p-( 1)  4p-( 1)
               1/2      1/2
                           0+
  1s ( 2)  2s ( 1)  4s ( 1)
               1/2      1/2
                           0+
 *
  1s ( 2)  2p ( 1)  4d-( 1)
               3/2      3/2
                           0-
  1s ( 2)  2s ( 1)  4p-( 1)
               1/2      1/2
                           0-
  1s ( 2)  2p-( 1)  4s ( 1)
               1/2      1/2
                           0-
```
where `CSF(s)` is in *jj*-notation, which corresponds to `name` entry of `.en` file.

Before running `jj2LSJ` module in GRASP2K, we need to correct some bugs in their source code.


# Fix GRASP2K source code and Install
## Download

The source code can be found in http://cpc.cs.qub.ac.uk/summaries/ADZL_v1_1.html

Extract it in an appropriate location.

## Bug fix
There are (at least) two bugs to run with FAC output.

1. In the first few lines of `asf2ls` subroutine of the `jj2lsj_code.f90`, shown below. The bold are what need to be added to the index of `ISPAR` array.
<pre><code>
wb = zero
do LS_number = 1, asf_set_LS%csf_set_LS%nocsf
   if ((asf_set_LS%csf_set_LS%csf(LS_number)%parity == "+" &
      .and.  ISPAR(iw1+<b>NCFMIN-1</b>) == 1)  .or.                        &
      (asf_set_LS%csf_set_LS%csf(LS_number)%parity  == "-" &
      .and.  ISPAR(iw1+<b>NCFMIN-1</b>) == -1)) then               
</code></pre>

2. In `setLS_job_count()` subroutine of the `jj2lsj_code.f90`, the following block inside the `if (all_occupation(isubc).eq.0) then` needs to be changed to the below.
<pre><code>
if(all_occupation(isubc).eq.0) then         
   if(isubc.gt.1) then
      Li(isubc)  = 0;          L_i(isubc) = L_i(isubc-1)
      Si(isubc)  = 0;          S_i(isubc) = S_i(isubc-1)
   else
      Li(isubc)  = 0;          L_i(isubc) = 0
      Si(isubc)  = 0;          S_i(isubc) = 0
   endif
   if(isubc .lt. asf_set_LS%csf_set_LS%nwshells) then
      call setLS_job_count(isubc + 1, rez)
   else
      if(ittk(S_i(isubc),L_i(isubc),J).eq.1)          &
           call setLS_action(action_type, rez) !rez=rez+1
   end if
else
</code></pre>


## Compilation and installation

Go to root directory of `GRASP2K`, run
```
source make-environment-gfort
```
(for gfortran compiler. If you are using another compiler, see README of `GRASP2K`.)

Go to `src` directory and run `make`
```
cd src
make
```
It will generate `jj2lsj` executable in `$GRASP$/bin` directory.


# Run jj2LSJ
## Preparation
You need `isodata` file in the current directory to run `jj2lsj`.
An example of `isodata` looks like this

```
Atomic number:
   4.0000000000000000     
Mass number (integer) :
   9.0000000000000000     
Fermi distribution parameter a:
  0.52338755531043146     
Fermi distribution parameter c:
   1.6101584296199354     
Mass of nucleus (in amu):
   1.5000000000000000     
Nuclear spin (I) (in units of h / 2 pi):
   0.0000000000000000     
Nuclear dipole moment (in nuclear magnetons):
   0.0000000000000000     
Nuclear quadrupole moment (in barns):
   0.0000000000000000     
```

I don't think the value itself is important, but `jj2lsj` checks its existence before running.

## Run
```
====================================================
       jj2lsj: Execution Begins ...
====================================================

  jj2lsj: Transformation of ASFs from a jj-coupled CSF basis
          into an LS-coupled CSF basis  (Fortran 95 version)
          (C) Copyright by   G. Gaigalas and Ch. F. Fischer,
          (2011).

Name of state
```
If you run `jj2lsj`, the above header will be printed.
Since `jj2lsj` is an interactive software, we need to input each line through shell.

First, we need tell the path to the `*.c` and `*.cm` files (they should be in the same directory).
For the second block in the above example, input `Be_002` in the shell.

```
Be_002
 Loading Configuration Symmetry List File ...
 There are 9 relativistic subshells;
 There are 36 relativistic CSFs;
  ... load complete;

 Mixing coefficients from a CI calc.?
```
Then it will ask us some questions. Probably we would simply hit `y` several times.

Finally, if it completes, something like below will be printed.
```
====================================================
       jj2lsj: Execution Finished ...
====================================================
Wall time:
     515 seconds

Finish Date and Time:
  Date (Yr/Mon/Day): 2018/11/12
  Time (Hr/Min/Sec): 08/06/13.467
  Zone: +0100

jj2lsj: Execution complete.
```

Then you get a file `Be_002.lsj.lbl`, which shows compositions of LSJ basis for each level.
The first several lines of the output will look like
```
Pos   J   Parity      Energy Total      Comp. of ASF
47    0     +           -28.521378592     100.000%
        0.99710485    0.99421808   1s(2).2s_2S.4s_1S
        0.07603890    0.00578191   1s(2).2p_2P.4p_1S
69    0     +           -28.367273893     100.000%
       -0.99996380    0.99992759   1s(2).2p_2P.4p_3P
97    0     +           -28.358234660     100.000%
        0.99706878    0.99414616   1s(2).2p_2P.4p_1S
       -0.07603578    0.00578144   1s(2).2s_2S.4s_1S

```

The first line shows the composition of the first CSF basis of `*.c` file,
```
1s ( 2)  2p ( 1)  4p ( 1)
             3/2      3/2
                         0+
```
The second-third line shows its composition.
In this case, this state consists of 99.7 \% of 2s4s 1S and 0.6 \% of 2p4p 1S.
