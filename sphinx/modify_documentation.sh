

rm -fR _build/*
rm -fR ../docs/*
#sphinx-apidoc -f -o autodoc/ ../k_quant/
make html
cp -va _build/html/* ../docs/
touch ../docs/.nojekyll
git add ..
