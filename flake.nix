{
  inputs.nixpkgs.url = "github:NixOS/nixpkgs/nixos-unstable";
  inputs.py-utils.url = "github:rmcgibbo/python-flake-utils/2arg";
  inputs.utils.url = "github:numtide/flake-utils";

  outputs = { self, nixpkgs, utils, py-utils }: rec {
    overlay = py-utils.lib.mkPythonOverlay (pySelf: pySuper: {
      debugpy = pySuper.debugpy.overridePythonAttrs (old: {
        doCheck = !pySelf.python.isPy38;
      });
      pyzmq = pySuper.pyzmq.overridePythonAttrs (old: {
        doCheck = !pySelf.python.isPy38;
      });
      mdtraj = pySelf.callPackage ./. {
        # right now, only build docs on python39
        buildDocs = !pySelf.python.isPy310;
      };
      mdtraj-generic = pySelf.callPackage ./. {
        buildDocs = false;
        enableNativeVectorIntrinsics = false;
      };
    });
  } //
  utils.lib.eachDefaultSystem (system:
    let
      pkgs = import nixpkgs {
        inherit system;
        overlays = [ self.overlay ];
      };
    in
    {
      packages = pkgs;
    });
}
