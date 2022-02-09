{
  inputs.nixpkgs.url = "github:NixOS/nixpkgs/nixos-unstable";
  inputs.py-utils.url = "github:rmcgibbo/python-flake-utils";
  inputs.utils.url = "github:numtide/flake-utils";

  outputs = { self, nixpkgs, utils, py-utils }: {
    overlay = py-utils.lib.mkPythonOverlay (pkgs: {
      mdtraj = pkgs.callPackage ./. {
        # right now, only build docs on python39
        buildDocs = !pkgs.python.isPy310;
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
