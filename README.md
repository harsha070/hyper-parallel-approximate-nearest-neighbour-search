# Hyper-Parallel Approximate Nearest Neighbor Search

### Project
Welcome to our project! We have implemented a parallelized version of approximate nearest neighbor search (ANN). Our code builds trees across different nodes and utilized multi-threading within each node. We have tested the implementation on Harvard's academic cluster. An analysis of the speedup can be found in the `report.pdf` file. This file describes our algorithm in depth, analyses the achieved speedup and gives an outlook on future work. We hope that you enjoy our project!

### Structure
The main function to create a search tree is located in `src/core/SearchTree.cpp`. The code to build a forest of trees is located in `src/core/SearchTree.cpp`. Header files for the SearchTree and the SearchForest are in the same folder. If you want to track the progress of our project have a look into the milestone folder. Here we share the presentations of different stages of building our application. The `scaling` folder contains files to recreate the weak and strong scaling analysis found in the report. Note that running the `./src/test.cpp file` just utilizes that available threads on the current machine and is intended to check correctness of the code. To evaluate performance, please run the strong and weak scaling job scripts in the `scaling` folder.

### Compilation
To test the tree building and inference run the test file in `src`. This will build a tree from random points and perform inference for a point to check the correct construction of the index. 

#### Building Search Tree

On the academic cluster, run the following commands

```bash
cd src/
spack load gcc
spack load openmpi@4.1.6
make test
```

#### Building Search Forest
On the academic cluster, run the following commands

```bash
cd src/
spack load gcc
spack load openmpi@4.1.6
make main
```

### Runtime

#### Building Search Tree

```bash
./test
```

#### Building Search Forest

```bash
./main
```
