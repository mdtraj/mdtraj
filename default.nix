{ buildPythonPackage
, fetchPypi
, setuptools
, cython
, numpy
, zlib
, pyparsing
, astunparse
, scipy
, networkx
, matplotlib
, scikitlearn
, pandas
, tables
, psutil
, gsd
, numpydoc
, nbformat
, nbconvert
, sphinx_rtd_theme
, sphinx
, pytestCheckHook
, pytest-xdist
,}:

let
  msmb-theme = buildPythonPackage rec {
    pname = "msmb_theme";
    version = "1.2.0";
    src = fetchPypi {
      inherit pname version;
      sha256 = "0b77yjk5q8kdp7bdlqlwi33hjirzp7bbblpqi4a2gy1an3ijzp4v";
    };
    propagatedBuildInputs = [
      numpydoc
      sphinx
      sphinx_rtd_theme
    ];
    # no tests
    doCheck = false;    
  };
in
buildPythonPackage {
  pname = "mdtraj";
  version = "0.0";
  src = ./.;

  buildInputs = [
    setuptools
    cython
    numpy
    zlib
  ];

  propagatedBuildInputs = [
    pyparsing
    astunparse
    scipy
    networkx
  ];

  checkInputs = [
    pytest-xdist
    pytestCheckHook

    matplotlib
    scikitlearn
    pandas
    tables
    psutil
    gsd

    # for docs
    nbformat
    nbconvert
    sphinx
    msmb-theme
  ];

  disabledTests = [
    # These tests require network access
    "test_pdb_from_url"
    "test_1vii_url_and_gz"
  ];

  # Ensure mdconvert is on the PATH and don't
  # let us import from the src directory
  preCheck = ''
    export PATH=$out/bin:$PATH
    rm -rf mdtraj/
  '';

  postInstall = ''
    mkdir -p $out/share/docs
    (cd docs && make html && cp -r _build/html $out/share/docs)
  '';
}
