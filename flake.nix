{
  description = "MFI — Modern Fortran Interfaces (BLAS/LAPACK)";

  inputs = {
    # Pin to 24.11 — same as the old shell.nix files
    # This gives us CUDA 12.3 (gpu-modern) and 11.8 (gpu-legacy)
    nixpkgs.url = "github:NixOS/nixpkgs/nixos-24.11";
    flake-utils.url = "github:numtide/flake-utils";
  };

  outputs = { self, nixpkgs, flake-utils }:
    flake-utils.lib.eachDefaultSystem (system:
      let
        # Overlay: fortran-fpm 0.13.0 (PR #506818)
        # Remove this once merged into nixpkgs
        fpmOverlay = final: prev: {
          fortran-fpm = prev.fortran-fpm.overrideAttrs (old: rec {
            version = "0.13.0";
            src = prev.fetchurl {
              url = "https://github.com/fortran-lang/fpm/releases/download/v${version}/fpm-${version}.F90";
              hash = "sha256-ABz/bPEUXyFbqgiIuieswGzqMKibedGovpfbP/+8jNI=";
            };
          });
        };

        pkgs = import nixpkgs {
          inherit system;
          config.allowUnfree = true;
          overlays = [ fpmOverlay ];
        };

        cpuLibs = [
          pkgs.hdf5
          pkgs.hdf5-fortran
          pkgs.blas
          pkgs.lapack
        ];

        cudaModern = pkgs.cudaPackages_12_3;
        cudaModernLibs = [
          cudaModern.libcublas
          cudaModern.libcublas.dev
          cudaModern.cuda_cudart
          cudaModern.cuda_cudart.dev
          cudaModern.cuda_nvcc
          cudaModern.cuda_cccl
        ];

        # Legacy: only CUDA/driver is old (11.8), rest is current nixpkgs
        cudaLegacy = pkgs.cudaPackages_11_8;
        cudaLegacyLibs = [
          cudaLegacy.libcublas
          cudaLegacy.libcublas.dev
          cudaLegacy.libcusolver
          cudaLegacy.libcusolver.dev
          cudaLegacy.cuda_cudart
          cudaLegacy.cuda_cudart.dev
          cudaLegacy.cuda_nvcc
        ];

        # Wrapper: provide `fpm` command pointing to fortran-fpm
        fpmAlias = pkgs.writeShellScriptBin "fpm" ''
          exec ${pkgs.fortran-fpm}/bin/fortran-fpm "$@"
        '';

        commonBuildInputs = [
          pkgs.pkg-config
          pkgs.gfortran
          pkgs.fortran-fpm
          pkgs.fypp
          fpmAlias
        ];

        mkGpuShell = { cudaLibs }: pkgs.mkShell {
          nativeBuildInputs = commonBuildInputs;
          buildInputs = cpuLibs ++ cudaLibs;
          shellHook = ''
            export CPATH="${pkgs.lib.makeSearchPath "include" cudaLibs}:$CPATH"
            export LIBRARY_PATH="${pkgs.lib.makeLibraryPath (cpuLibs ++ cudaLibs)}:$LIBRARY_PATH"
            export LD_LIBRARY_PATH="${pkgs.lib.makeLibraryPath (cpuLibs ++ cudaLibs)}:$LD_LIBRARY_PATH"
          '';
        };

        mkCpuShell = pkgs.mkShell {
          nativeBuildInputs = commonBuildInputs;
          buildInputs = cpuLibs;
          shellHook = ''
            export LIBRARY_PATH="${pkgs.lib.makeLibraryPath cpuLibs}:$LIBRARY_PATH"
            export LD_LIBRARY_PATH="${pkgs.lib.makeLibraryPath cpuLibs}:$LD_LIBRARY_PATH"
          '';
        };
      in
      {
        devShells = {
          cpu-only = mkCpuShell;
          gpu-modern = mkGpuShell { cudaLibs = cudaModernLibs; };
          gpu-legacy = mkGpuShell { cudaLibs = cudaLegacyLibs; };
          default = mkCpuShell;
        };
      }
    );
}
