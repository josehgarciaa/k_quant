

rm autodoc/*
sphinx-apidoc -o autodoc/ ../kquant
make html
cp -a _build/html ../docs/
git add ..
