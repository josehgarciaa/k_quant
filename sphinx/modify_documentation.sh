

rm -f autodoc/*
sphinx-apidoc -o autodoc/ ../k_quant
make html
cp -a _build/html ../docs/
git add ..