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

        # ZLUDA: drop-in CUDA/cuBLAS replacement for AMD GPUs.
        # ZLUDA is a pre-built binary (not compiled from nixpkgs — pkgs.zluda pulls in
        # rocmlir-rock which is broken in nixpkgs 24.11).  Download it from:
        #   https://github.com/vosen/ZLUDA/releases
        # and point ZLUDA_PATH at its directory before entering this shell.
        #
        # Everything else IS provided by Nix:
        #   rocmPackages.clr         — HIP runtime: libamdhip64.so (ZLUDA needs this at run time)
        #   rocmPackages.rocm-runtime — HSA runtime: libhsa-runtime64.so
        #   cudaModern.{libcublas,cuda_cudart}.dev — CUDA headers needed at compile time
        #
        # The only remaining host requirement is the AMD GPU kernel driver
        # (the amdgpu kernel module + firmware), which Nix cannot deliver.
        rocmLibs = [
          pkgs.rocmPackages.clr           # HIP runtime: libamdhip64.so + HIP headers
          pkgs.rocmPackages.rocm-runtime  # HSA runtime: libhsa-runtime64.so
        ];
        mkZludaShell = pkgs.mkShell {
          nativeBuildInputs = commonBuildInputs;
          # CPU libs + CUDA headers (compile time) + ROCm/HIP stack (ZLUDA runtime deps)
          buildInputs = cpuLibs ++ rocmLibs ++ [
            cudaModern.libcublas.dev
            cudaModern.cuda_cudart.dev
            cudaModern.cuda_cccl
          ];
          shellHook = ''
            # CUDA headers so cuda_runtime.h / cublas_v2.h are found at compile time
            export CPATH="${pkgs.lib.makeSearchPath "include" [
              cudaModern.libcublas.dev
              cudaModern.cuda_cudart.dev
              cudaModern.cuda_cccl
            ]}:$CPATH"
            # ROCm/HIP stack so ZLUDA can resolve libamdhip64/libhsa-runtime64 at run time
            export LIBRARY_PATH="${pkgs.lib.makeLibraryPath (rocmLibs ++ cpuLibs)}:$LIBRARY_PATH"
            export LD_LIBRARY_PATH="${pkgs.lib.makeLibraryPath (rocmLibs ++ cpuLibs)}:$LD_LIBRARY_PATH"
            # Wire in the user-supplied ZLUDA directory (must come first to override stubs)
            if [ -n "$ZLUDA_PATH" ]; then
              export LIBRARY_PATH="$ZLUDA_PATH:$LIBRARY_PATH"
              export LD_LIBRARY_PATH="$ZLUDA_PATH:$LD_LIBRARY_PATH"
              echo "ZLUDA shell ready (ZLUDA_PATH=$ZLUDA_PATH)."
              echo "  Build: make && fpm build --profile zluda"
              echo "  Run:   MFI_USE_CUBLAS=1 ./build/gfortran_*/app/app"
            else
              echo "WARNING: ZLUDA_PATH is not set."
              echo "  Download ZLUDA from https://github.com/vosen/ZLUDA/releases"
              echo "  then re-enter with: ZLUDA_PATH=/path/to/zluda nix develop .#gpu-zluda"
            fi
          '';
        };
      in
      {
        devShells = {
          cpu-only = mkCpuShell;
          gpu-modern = mkGpuShell { cudaLibs = cudaModernLibs; };
          gpu-legacy = mkGpuShell { cudaLibs = cudaLegacyLibs; };
          gpu-zluda  = mkZludaShell;
          default = mkCpuShell;
        };
      }
    );
}
