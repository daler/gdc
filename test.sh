cwd=`pwd`
cd doc && make doctest html &&
cd $cwd
nosetests --with-doctest --doctest-extension=.rst
