

rm -fR build/*
rm -fR ../docs/*
#sphinx-apidoc -f -o autodoc/ ../k_quant/
make clean
make html
cp -va build/html/* ../docs/
touch ../docs/.nojekyll
git add ..
