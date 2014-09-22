set -x

download_and_unpack() {
    local url=${1}
    local md5=${2}
    local unpack=${3}
    local file=`basename $url`
    
    wget $url
    local received_md5=`md5sum "$file" | awk '{print $1}'`

    if [ "$md5" != "$received_md5" ] ; then
        echo 'MD5 mismatch!'
        exit 1
    fi
    
    tar -xzvf $file
}

# Fetch ASTOR 0.4
url="https://pypi.python.org/packages/source/a/astor/astor-0.4.tar.gz"
reference_md5='32ee1e50e88cccfc56315b19b09449d3'
unpack="astor-0.4"
target="astor"
download_and_unpack $url $reference_md5 $unpack
rm -rf $target
mkdir $target
cp -r $unpack/astor/* $target/
cp $unpack/LICENSE $target/
rm -rf $unpack
rm -rf `basename $url`

# Fetch PyParsing 2.0.2
url="https://pypi.python.org/packages/source/p/pyparsing/pyparsing-2.0.2.tar.gz"
reference_md5="b170c5d153d190df1a536988d88e95c1"
unpack="pyparsing-2.0.2"
target="pyparsing.py"
download_and_unpack $url $reference_md5 $unpack
cp $unpack/pyparsing.py $target
rm -rf $unpack
rm -rf `basename $url`