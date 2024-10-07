# Topic
As discussed in [issue #1](https://github.com/jooooow/rtfmm/issues/1), regularization by copying image cells will introduce inevitable error.

Therefore, we try to combine:
+ FMM
+ Ewald
+ Regularization

to sovle the problem.

# Build
```shell
bash ./script/build.sh
```

# Execution

+ test_fmm : kernel-independent FMM without PBC
+ test_pbc : kernel-independent FMM with PBC
+ test_reg : kernel-independent FMM with PBC and regularization
    - WARNING : currently, enable PBC and Regularization at the same will increase the error
+ test_ef : kernel-independent FMM with ewald-based-PBC and regularization


# Usage

### syntax
```shell
python3 script/run.py [exec] <options>
```
### example
```shell
python3 script/run.py test_fmm -n 24000 -a fd -P 8 -i 0 -r 0 --seed=1234 -v 1 
```

# Output

+ `./result` will be automatically created after running the script.

+ `./result/[time]@[exec].log` is the log file of the execution.
