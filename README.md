# LINE: Loop Integrals Numerical Evaluation

**LINE** (which stands for Loop Integrals Numerical Evaluation), is a
tool to compute **Feynman integrals** by numerically solving
differential equations via series expansion. **LINE** is written in `C`
/ `C++` and leverages the well-known **GMP** family of libraries for
**arbitrary precision arithmetic**, aiming to achieve efficiency and
accessibility in order to go beyond proof of concept and make
large-scale cluster computations more feasible.

## Dependencies

**LINE** depends on the following libraries:

-   **GMP (GNU Multiple Precision Arithmetic Library)**: 6.3.0
-   **MPFR (Multiple Precision Floating-Point Reliable Library)**: 4.2.1
-   **MPC (Multiple Precision Complex Library)**: 1.3.1
-   **OpenSSL**: 3.4.0
-   **Kira**: 2.3 (optional);

While other versions may work, the code has been tested and validated
only with the versions listed above. Using different versions may lead
to unexpected behavior or compatibility issues.

Kira is only needed if the user wants to generate boundaries using the
**LINE** implementation of the
[AMFlow](https://inspirehep.net/literature/1639025) method.

### Install GMP

To install GMP, follow these steps:

1.  Download the source code from the [official GMP web
    site](https://gmplib.org) (e.g.
    [gmp-6.3.0.tar.xz](https://gmplib.org/download/gmp/gmp-6.3.0.tar.xz)).
2.  Navigate to the directory containing the source code, then execute:

``` bash
./configure
make
make check
sudo make install
```

### Install MPFR

MPFR depends on GMP, so ensure that GMP is installed first. Then:

1.  Download the source code from the [official MPFR web
    site](https://www.mpfr.org) (e.g.
    [mpfr-4.2.1.tar.xz](https://www.mpfr.org/mpfr-current/mpfr-4.2.1.tar.xz)).
2.  Proceed as for the MPFR installation. If GMP is installed in a
    non-standard directory, or if more than one version is present,
    specify the paths to GMP using `--with-gmp-include` and
    `--with-gmp-lib` to indicate the folders containing the header file
    `gmp.h` and the library `libgmp.so`, respectively (e.g.
    `/usr/local/include` and `/usr/local/lib`, respectively):

``` bash
./configure --with-gmp-include=/usr/local/include --with-gmp-lib=/usr/local/lib
make
make check
sudo make install
```

### Install MPC

MPC requires both GMP and MPFR. To install it:

1.  Download the source code from the [official MPC web
    site](https://www.multiprecision.org) (e.g.
    [mpc-1.3.1.tar.gz](https://ftp.gnu.org/gnu/mpc/mpc-1.3.1.tar.gz)).
2.  Proceed as for the MPC installation, using `--with-mpfr-include` and
    `--with-mpfr-lib~` to specify include and lib directories if
    necessary:

``` bash
./configure --with-gmp-include=/usr/local/include --with-gmp-lib=/usr/local/lib --with-mpfr-include=/usr/local/include --with-mpfr-lib=/usr/local/lib
make
make check
sudo make install
```

### Install OpenSSL

The required libraries are `libcrypto` and `libssl`, which are part of
OpenSSL.

On **Linux**, OpenSSL can be installed with:

``` bash
sudo apt-get install libssl-dev
```

On **MacOs**, you can install OpenSSL with **Homebrew**:

``` bash
brew install openssl
```

### Install Kira

Kira can be installed following the instructions on the [project
official git repository](https://gitlab.com/kira-pyred/kira). To run
Kira, [Fermat](http://home.bway.net/lewis/) is required. Make sure to
set the environment variable `FERMATPATH` to indicate the path to the
Fermat executable, for instance adding:

``` bash
export FERMATPATH="/path/to/fer64"
```

to your shell\'s initialization file (e.g., `.bashrc`, `.bash_profile`,
etc.).

## Compile and check

The source code can be compiled with `make`. A suit of light tests is
provided in the `check/` folder to check that everything has been
properly installed. The tests can be run through the script
`check_run.sh`:

``` bash
make
./chek_run.sh
```

To perform multiple checks simultaneously using N CPU cores, run
`check_run.sh` with the `--nthreads N` option (or `-n N` for short):

``` bash
./chek_run.sh -n N
```

In order to check the `Kira` installation as well (together with the
standard tests), run with the `--kira` option (or `-k` for short):

``` bash
./chek_run.sh -k
```

### More tests

Further tests will be made available soon. At present, these files have
not been included in this repository because the DE matrices and the
validation files used to verify the correctness of intermediate
calculation results require significant storage space. To address this,
a separate Git repository will be created to host these files. Users
will be able to download the repository and integrate it as a Git
submodule into their local setup. Stay tuned for updates on this
resource.

## Usage

**LINE** requires the following input files:

-   The Differential Equation (DE) matrices.
-   The boundary values for the Master Integrals (MIs) to be propagated.
-   Info about the branch cuts for the considered topology.
-   An input card indicating the target point, the starting point, the
    required number of digits for the result as well as the number of
    orders for its Laurent expansion in the dimensional regularization
    parameter epsilon.

If you want to use the **LINE** implementation of the AMFlow method, you
must also provide:

-   A list of master integrals.
-   Information about the topology.

Note that if you do not have the DEs with respect to the kinematic
invariants, the AMFlow method can still be used by providing the
integrals you want to calculate in the list of MIs.

### Structure

Each topology is associated with a working directory `work-dir/`. The
topology-related files (everything except the input card) must be stored
into the `common/` subdirectory. After the first run, the subdirectory
`points/` is also created. Each time a run produces results for a new
phase-space point, a directory is created within `points/`, where the
calculated values of the MIs are saved. When propagating to a new point,
**LINE** automatically checks if the starting point has been used
previously as a target. If so, it can reuse the result from the previous
run as the boundary for the new propagation.

The typical structure of the working directories is as follows:

    work-dir/
    ├── common/  # run-independent files
    │   ├── 0.txt  # DE matrix w.r.t. 1st invariant
    │   ├── 1.txt  # DE matrix w.r.t. 2nd invariant
    │   ├── 2.txt  # DE matrix w.r.t. 3rd invariant
    │   ├── ...
    │   ├── vars.txt  # list of symbols for the invariants
    │   ├── MIs.txt  # list of MIs
    │   ├── topology.txt  # topology info
    │   └── branch_cuts.txt  # list of branch cuts invariants
    └── points/
        ├── hash-target/  # some target point
        │   ├── sol/  # solution folder
        │   ├── hash-start0/  # some starting point
        │   │   └── ...  # cache files for starting-target propagation
        │   ├── hash-start1/  # another starting point
        │   │   └── ...
        │   └── ...
        ├── hash-target1/  # another target point
        │   ├── sol/
        │   ├── hash-start0/
        │   │   └── ...
        │   ├── hash-start3/
        │   │   └── ...
        │   └── ...
        └── ...

### How to run

A typical example of an input card looks as follows:

    work-dir: 1L-triangle-full  # working-directory
    loops: 1  # number of loops
    order: 6  # number of order for the Laurent series beyond the leading coefficient
    precision: 16  # number of digits for the result
    point: [  # target point
      s = 1,
      p12 = 2,
      p22 = -1/3,
      m12 = (1 -1),  # (real-part imaginary-part)
      m22 = (8/3 -2),
      m32 = (17 -1/4)
    ]
    starting-point: [  # starting point
      s = 50,
      p12 = 2,
      p22 = -1/3,
      m12 = 5,
      m22 = 7,
      m32 = 10
    ]

To use the **LINE** implementation of the AMFlow method (only
real-valued kinematics), there is no need to specify a starting point,
while the options `exit sing: -1` and `gen-bound 1` must be used:

    work-dir: 1L-triangle-full
    loops: 1
    order: 6
    precision: 16
    point: [
      s = 50,
      p12 = 2,
      p22 = -1/3,
      m12 = 5,
      m22 = 7,
      m32 = 10
    ]
    exit-sing: -1
    gen-bound: 1

You can run with:

``` bash
./line -i path/to/input/card -r path/to/result/file > out.log
```

To write the computed results to the cache folders (so that they can be
used as a boundary for later propagations) use `--write 1` (`-w 1` for
short). If instead, for a given target point, you want to check that the
result is in agreement with a previous run that used the same target but
a different starting point, use `--write 0`.

To use input boundary values for the MIs, use
`--bound-dir path/to/bound/dir`, where `path/to/bound/dir` is a folder
containing files with the values of the MIs in several epsilon values:

-   `bound0.txt`: values of the MIs computed in the 1st epsilon values;
-   `bound1.txt`: values of the MIs computed in the 2nd epsilon values;
-   ...

Each boundary file must be in the MPC format
`(real-part imaginary-part)`. For instance, for a problem with 7 MIs:

``` bash
(3.01388546321267...521e3 0)  # boundary value for the 1st MI computed in 1st epsilon
(4.21708961212546...355e3 0)  # boundary value for the 2nd MI computed in 1st epsilon
(6.01659689683808...329e2 0)  # boundary value for the 3rd MI computed in 1st epsilon
(6.02395241231279...583e2 2.00628403154449...309e0)  # ...
(6.01249274999828...526e2 0)
(-7.55253164143286...510e-2 -1.02172397265745...499e-1)
```

(digits are truncated with `...` just for ease of visualization). To
determine the list of required epsilon values (which depend on the
inputs `order` and `precision`), you can use the `--eps-list` option to
make `line` stop immediately after printing the epsilon list.

For some of the provided examples, automated boundaries can be generated
through **Expansion By Regions** (EBR) using `exit-sing: 1` in the input
card. When doing so, if no `starting-point` is provided in the input
card, the one specified in the file `work-dir/common/initial_point.txt`
is used. The other relevant files in the `common/` directory to generate
boundaries with EBR are:

-   `bound_behav.txt`: specify power behaviors of the MIs around the
    expansion point
-   `bound_build.txt`: instruction to build leading coefficients of MIs
    from simpler subtopologies (leading coefficients for the others are
    automatically obtained when solving the DE)

Please refer to the provided examples for the format of the input files
located in the `work-dir/common/` directory.

## Alpha version and future improvements

This code is an alpha version of the tool. While functional, there is
still significant room for improvement, particularly in terms of
efficiency, code compactness, and modularity. Many enhancements have
already been conceptualized and a beta version including these
improvements is expected to be released in the near future. Stay tuned
for updates!

## Reporting issues

Please note that as this is an alpha version, the code may be prone to
bugs. If you encounter any issues, we encourage you to [contact
us](mailto:renatomaria.prisco@unina.it) for support. Please include a
detailed description of the problem and steps to reproduce it. This will
help us to diagnose and resolve the issue more efficiently, improving
the tool for future releases.

## Citation

If you use **LINE** in your work, please cite our paper:

------------------------------------------------------------------------

R. M. Prisco, J. Ronca, F. Tramontano, *LINE: Loop Integrals Numerical
Evaluation* - [arXiv:2501.01943](https://arxiv.org/abs/2501.01943)

------------------------------------------------------------------------

This contributes to the continued development and improvement of the
project.
