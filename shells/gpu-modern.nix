{ pkgs ? import (fetchTarball "https://github.com/NixOS/nixpkgs/archive/refs/tags/24.11.tar.gz") { config.allowUnfree = true; } }:

let
  cudaPkgs = pkgs.cudaPackages_12_3;
  cudaLibs = [
    cudaPkgs.libcublas
    cudaPkgs.libcublas.dev
    cudaPkgs.cuda_cudart
    cudaPkgs.cuda_cudart.dev
    cudaPkgs.cuda_nvcc
    cudaPkgs.cuda_cccl
  ];
  allLibs = [
    pkgs.hdf5
    pkgs.hdf5-fortran
    pkgs.blas
    pkgs.lapack
  ] ++ cudaLibs;
in
pkgs.mkShell {
  nativeBuildInputs = [
    pkgs.pkg-config 
    pkgs.gfortran
  ];

  buildInputs = allLibs;

  shellHook = ''
    export CPATH="${pkgs.lib.makeSearchPath "include" cudaLibs}:$CPATH"
    export LIBRARY_PATH="${pkgs.lib.makeLibraryPath allLibs}:$LIBRARY_PATH"
    export LD_LIBRARY_PATH="${pkgs.lib.makeLibraryPath allLibs}:$LD_LIBRARY_PATH"
  '';
}
