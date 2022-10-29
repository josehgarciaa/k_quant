

sphinx-apidoc -f -o autodoc/ ../k_quant/
make html
cp -va _build/html/* ../docs/
git add ..
