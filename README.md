## Description

## Running without nix

First, make sure the submodules are available

```
git submodule update --init --recursive
```

The project can then be configured as you would with Kokkos. For example, with OpenMP

```
cmake -DCMAKE_BUILD_TYPE=Release -DKokkos_ENABLE_OPENMP=ON -B build
```

```
cmake --build build
```

```
./build/main
```

## Running with nix
