

rm -fR _build/*
rm -fR ../docs/*
touch _build/html/.nojekyll
sphinx-apidoc -f -o autodoc/ ../k_quant/
make html
cp -va _build/html/* ../docs/
git add ..
