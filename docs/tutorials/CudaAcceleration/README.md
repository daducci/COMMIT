# Enable CUDA GPU Acceleration

This tutorial illustrates how to enable the CUDA GPU acceleration for faster model fitting.

## Getting Started

This tutorial takes as starting point the [getting started](https://github.com/daducci/COMMIT/tree/master/docs/tutorials/GettingStarted) tutorial.

## Build the linear operator A on GPU

Once the getting started tutorial runs properly, focus on the commands `mit.set_threads()` and `mit.build_operator()`. The command `mit.set_threads()` allows to indicate (through the parameter `nthreads`) how many CPU threads will use COMMIT to build linear operator **A**. Then, the command `mit.build_operator()` builds **A** by using the same number of threads especified by `nthreads`. The GPU acceleration can be enabled by setting the parameter `nthreads` equal to zero in `mit.set_threads()`.

```python
nthreads = 0
mit.set_threads( nthreads=nthreads )
mit.build_operator()
```

Once these commands area executed, COMMIT will build the linear operator ***A** but now on the GPU. The output should be something similar to this:

```
-> Building linear operator A:
        * checking availability of CUDA...  [ OK ]
        * number of CUDA GPUs detected: 	1
        * using GPU with ID 0... 			[ GeForce RTX 2080 SUPER ]
        * using 0.31 GB of total 8.37 GB... [ OK ]
        * compute capability: 7.5 			[ OK ]
        * constant values ... 				[ OK ]
        * vectors x&y ... 					[ OK ]
        * pre-processing ... 				[ OK ]
        * loading LUTs ... 					[ OK ]
        * A  operator... 					[ OK ]
        * A' operator... 					[ OK ]
   [ 0.5 seconds ]
```

If there are more than one CUDA capable GPU installed, COMMIT will use by default the GPU with ID 0. The selected GPU can be changed with the parameter `select_gpu` in `mit.set_threads()`. For example, assuming there are at least two CUDA capable GPUs installed, the following commands build the linear operator **A** on the GPU with ID 1 rather than on the GPU with ID 0:

```python
nthreads = 0
gpu_id = 1
mit.set_threads( nthreads=nthreads, select_gpu=gpu_id )
mit.build_operator()
```

To show a list of GPUs and their IDs, open a system shell and run the command `nvidia-smi`. This command should output something similar to:

```
                             +-----------------------------------------------------------------------------+
                             | NVIDIA-SMI 375.82                 Driver Version: 375.82                    |
                             |-------------------------------+----------------------+----------------------+
                             | GPU  Name        Persistence-M| Bus-Id        Disp.A | Volatile Uncorr. ECC |
                             | Fan  Temp  Perf  Pwr:Usage/Cap|         Memory-Usage | GPU-Util  Compute M. |
                             |===============================+======================+======================|
this is the GPU ID  -------> |   0  GeForce GTX TIT...  Off  | 0000:05:00.0     Off |                  N/A |
                             | 22%   55C    P8    31W / 250W |  11853MiB / 12205MiB |      0%      Default |
                             +-------------------------------+----------------------+----------------------+
this is the GPU ID  -------> |   1  GeForce GTX TIT...  Off  | 0000:06:00.0     Off |                  N/A |
                             | 22%   60C    P8    18W / 250W |    114MiB / 12207MiB |      0%      Default |
                             +-------------------------------+----------------------+----------------------+
this is the GPU ID  -------> |   2  GeForce GTX TIT...  Off  | 0000:09:00.0     Off |                  N/A |
                             | 27%   66C    P2    72W / 250W |   8452MiB / 12207MiB |      0%      Default |
                             +-------------------------------+----------------------+----------------------+
```

**NOTE:** At this moment, COMMIT does not have support for multi-GPU acceleration.

## Clearing GPU memory

The commands apart from the command `mit.set_threads()` remain the same. But in the case when the GPU acceleration was enabled, the method `mit.A.destroy()` has to the executed in order to clear the GPU memory and reset the GPU for the next evaluation. That is, add the following command to the end of the script:

```python
if nthreads == 0:
	mit.A.destroy()
```

Then, something like the following is displayed:

```
-> Clearing GPU memory:
        * deleting A...   [ OK ]
        * deleting A'...  [ OK ]
        * deleting x&y... [ OK ]
        * deleting LUT... [ OK ]
        * reseting GPU... [ OK ]
```