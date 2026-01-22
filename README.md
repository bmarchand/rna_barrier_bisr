# rna_barrier_bisr_algorithms


This project presents Divide- and Subtree-algorithms, two techniques to solve the RNA energy barrier problem (that can be seen as the specialization of Bipartite Independent Set Reconfiguration under the TAR model) as presented in https://hal.science/hal-04094405/. The repository also contains the codes to produce the experiments and figures of the paper. All codes are in C++ and Python3 with the interface in Python3 only.


### Files and repositories

In order to avoid dependency management, we propose to use virtualenv as described below.
The project is entirely available on the current archive.

The different folders and files are:

	launcher_rna_bisr_solver.py: a parser to launch the tool in the command line;
    rna_barrier_bisr: the scripts for the two algorithms;
	experiments: the script to reproduce the experiments in https://hal.science/hal-04094405/;
	cpp: C++ routines;
    tests: a test suite to check the correctness of the program and setup.

### Launch

To set up the virtualenv environment in a shell, one can type:
```
virtualenv rna_barrier_bisr_env
source rna_barrier_bisr_env/bin/activate
pip install .
```
As a first example, one can put in the same shell: 
```bash
python3 launcher_rna_bisr_solver.py --method divide --s1 "..((((((..(.((((.((((........((((((...(((.......)))......)))))).......))))..)))))..(((((((((.(((.....))).)))))))))((((.((((((.....)))))).)))).((((((((....))))))))....))))))......((((((((((.....)))))))..)))...((((((((((.....))))))))))" --s2 ".....(((..(.((((.((((........((((((...(((.......)))......)))))).......))))..)))))..(((((((((.(((.....))).)))))))))((((.((((((.....)))))).)))).((((((((....))))))))....)))(((((((((((((((((((.....)))))))..))...))))))))))................"
```
* method is the chosen algorithm, either divide or subtree;
* s1 is the starting sequence and s2 the end sequence for the reconfiguration;
* The program should output the barrier for this instance in a few seconds.

The same instance can be solved with the other method by typing:
```bash
python3 launcher_rna_bisr_solver.py --method subtree --s1 "..((((((..(.((((.((((........((((((...(((.......)))......)))))).......))))..)))))..(((((((((.(((.....))).)))))))))((((.((((((.....)))))).)))).((((((((....))))))))....))))))......((((((((((.....)))))))..)))...((((((((((.....))))))))))" --s2 ".....(((..(.((((.((((........((((((...(((.......)))......)))))).......))))..)))))..(((((((((.(((.....))).)))))))))((((.((((((.....)))))).)))).((((((((....))))))))....)))(((((((((((((((((((.....)))))))..))...))))))))))................"
```

### Contributors

    Théo Boury
    Bertrand Marchand








