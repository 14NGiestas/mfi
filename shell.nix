# shell.nix
{ pkgs ? import <nixpkgs> {} }:
pkgs.mkShell {
  nativeBuildInputs = [
    pkgs.pkg-config 
    pkgs.fpm
  ];
  buildInputs = [
    pkgs.hdf5
    pkgs.hdf5-fortran
    pkgs.blas    # Provides blas.pc
    pkgs.lapack  # Provides lapack.pc
  ];
}

