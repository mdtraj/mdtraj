for i in `find . -name '*.py' -o -name '*.pyx'` # or whatever other pattern...
do
  if ! grep -q "Copyright" $i
  then
    cat LICENSE_HEADER $i >$i.new && mv $i.new $i
  fi
done