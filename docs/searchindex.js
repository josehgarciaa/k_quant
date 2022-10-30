Search.setIndex({"docnames": ["autodoc/k_quant", "autodoc/k_quant.visualization", "autodoc/modules", "example", "examples/example_0", "index", "modules"], "filenames": ["autodoc/k_quant.rst", "autodoc/k_quant.visualization.rst", "autodoc/modules.rst", "example.rst", "examples/example_0.rst", "index.rst", "modules.rst"], "titles": ["k_quant package", "k_quant.visualization package", "k_quant", "Tutorials", "A simple band structure calculation", "What is <span class=\"math notranslate nohighlight\">\\({\\boldsymbol k}\\)</span>\u00a0Quant?", "Documentation"], "terms": {"visual": [0, 2], "sciplot": [0, 2], "add_legend": [0, 1], "plot": [0, 1, 4], "plot4fig": [0, 1], "plot_zwitharrow": [0, 1], "class": [0, 4, 6], "bandstructur": [0, 2, 4, 5, 6], "lat_vec": [0, 4, 6], "h_k": [0, 6], "base": [0, 6], "object": [0, 6], "thi": [0, 4, 6], "comput": [0, 4, 5, 6], "differ": [0, 5, 6], "input": [0, 6], "model": [0, 4, 5, 6], "The": [0, 4, 6], "__init__": [0, 6], "method": [0, 6], "mai": [0, 6], "document": 0, "either": [0, 6], "level": [0, 5, 6], "docstr": [0, 6], "itself": [0, 6], "form": [0, 6], "i": [0, 4, 6], "accept": [0, 6], "two": [0, 6], "should": [0, 6], "mix": [0, 6], "choos": [0, 4, 6], "one": [0, 6], "convent": [0, 6], "consist": [0, 4, 6], "note": [0, 6], "do": [0, 6], "includ": [0, 6], "self": [0, 6], "paramet": [0, 6], "arg": [0, 1, 6], "section": [0, 6], "msg": [0, 6], "str": [0, 6], "human": [0, 6], "readabl": [0, 6], "string": [0, 6], "describ": [0, 6], "except": [0, 6], "code": [0, 6], "int": [0, 6], "option": [0, 1, 6], "error": [0, 6], "attribut": [0, 6], "momentum_rec2absmatrix": [0, 2, 6], "xlabel": [0, 2, 4, 6], "return": [0, 4, 6], "x": [0, 1, 6], "axi": [0, 6], "label": [0, 4, 6], "associ": [0, 6], "bandpath": [0, 2, 4, 6], "list": [0, 4, 6], "_description_": [0, 1, 6], "xaxi": [0, 2, 4, 6], "uniqu": [0, 6], "respect": [0, 6], "distanc": [0, 6], "between": [0, 4, 6], "k": [0, 4, 6], "point": [0, 4, 6], "ndarrai": [0, 6], "an": [0, 5, 6], "arrai": [0, 6], "band_kpoint": [0, 2, 6], "absolute_coord": [0, 6], "true": [0, 6], "defin": [0, 4, 6], "bool": [0, 6], "determin": [0, 6], "wherea": [0, 6], "ar": [0, 5, 6], "given": [0, 6], "cartesian": [0, 6], "unit": [0, 6], "2pi": [0, 6], "reciproc": [0, 6], "default": [0, 1, 6], "fals": [0, 6], "_summary_": [0, 1, 6], "_type_": [0, 1, 6], "compute_band": [0, 2, 4, 6], "fermi_energi": [0, 6], "0": [0, 4, 6], "proj_op": [0, 4, 6], "none": [0, 1, 6], "float": [0, 6], "ax": [0, 1, 4, 6], "plot_proj": [0, 6], "proj_rang": [0, 6], "hamiltonian_k": [0, 2, 6], "operator_k": [0, 2, 6], "oper": [0, 6], "kpoint": [0, 4, 6], "set_bandpath": [0, 2, 4, 6], "set_hamiltonian_k": [0, 2, 6], "toabsolutecoord": [0, 2, 6], "energies_inwindow": [0, 2, 6], "energy_window": [0, 6], "create_2dkgrid": [0, 2], "kwindow": 0, "npoint": [0, 1], "1": [0, 4], "legend": 1, "y": 1, "z": [1, 4], "kwarg": 1, "xx": 1, "yy": 1, "zz": 1, "uu": 1, "vv": 1, "zlabel": 1, "mask": 1, "3": [1, 4], "packag": [2, 5], "subpackag": 2, "submodul": 2, "modul": [2, 5], "content": 2, "band": [2, 3, 5], "densiti": [2, 5], "A": [3, 5], "simpl": [3, 5], "structur": [3, 5], "calcul": [3, 5], "For": 4, "exampl": [4, 5], "let": 4, "u": 4, "first": [4, 5], "hamiltonian": [4, 5], "momentum": [4, 5], "space": [4, 5], "which": [4, 5], "we": 4, "us": 4, "In": 4, "case": 4, "prototyp": 4, "p_z": 4, "electron": 4, "graphen": 4, "within": 4, "nearest": 4, "neighbor": 4, "approxim": 4, "lattic": 4, "vector": 4, "variabl": 4, "from": [4, 5], "numpi": 4, "import": 4, "sqrt": 4, "exp": 4, "dot": 4, "conj": 4, "min": 4, "max": 4, "ab": 4, "2": 4, "def": 4, "a_0": 4, "a_1": 4, "a2": 4, "hop": 4, "8": 4, "f_k": 4, "1j": 4, "Then": 4, "need": 4, "tupl": 4, "its": 4, "posit": 4, "normal": 4, "brilluoin": 4, "zone": 4, "number": 4, "next": 4, "onc": 4, "pass": 4, "k_quant": 4, "requir": 4, "function": [4, 5], "To": 4, "desir": 4, "path": 4, "you": 4, "npt": 4, "100": 4, "111": 4, "g": 4, "35": 4, "m": 4, "55": 4, "simmpl": 4, "sigma_x": 4, "sigma_i": 4, "final": 4, "result": 4, "matplotlib": 4, "pyplot": 4, "plt": 4, "fig": 4, "gcf": 4, "gca": 4, "proj_band": 4, "": [4, 5], "20": 4, "c": 4, "im": 4, "scatter": 4, "cmap": 4, "coolwarm": 4, "vmin": 4, "vmax": 4, "colorbar": 4, "set_xtick": 4, "set_xticklabel": 4, "savefig": 4, "proj_band_sigma_x": 4, "pdf": 4, "show": 4, "There": 5, "multipl": 5, "situat": 5, "where": 5, "research": 5, "creat": 5, "exploit": 5, "bloch": 5, "theorem": 5, "symmetri": 5, "intuit": 5, "howev": 5, "just": 5, "step": 5, "long": 5, "journei": 5, "toward": 5, "gener": 5, "quantit": 5, "predict": 5, "measur": 5, "quantiti": 5, "high": 5, "python": 5, "aim": 5, "provid": 5, "mani": 5, "spectral": 5, "properti": 5, "spin": 5, "orbit": 5, "project": 5, "conduct": 5, "hall": 5, "state": 5}, "objects": {"": [[0, 0, 0, "-", "k_quant"]], "k_quant": [[0, 0, 0, "-", "bands"], [0, 0, 0, "-", "densities"], [1, 0, 0, "-", "visualization"]], "k_quant.bands": [[0, 1, 1, "", "bandstructure"], [0, 3, 1, "", "energies_inwindow"]], "k_quant.bands.bandstructure": [[0, 2, 1, "", "Momentum_Rec2AbsMatrix"], [0, 2, 1, "", "XLabels"], [0, 2, 1, "", "Xaxis"], [0, 2, 1, "", "band_kpoints"], [0, 2, 1, "", "bandpath"], [0, 2, 1, "", "compute_bands"], [0, 2, 1, "", "hamiltonian_k"], [0, 2, 1, "", "operator_k"], [0, 2, 1, "", "set_bandpath"], [0, 2, 1, "", "set_hamiltonian_k"], [0, 2, 1, "", "toAbsoluteCoords"]], "k_quant.densities": [[0, 3, 1, "", "create_2Dkgrid"], [0, 3, 1, "", "energies_inwindow"]], "k_quant.visualization": [[1, 0, 0, "-", "sciplots"]], "k_quant.visualization.sciplots": [[1, 3, 1, "", "add_legends"], [1, 3, 1, "", "plot"], [1, 3, 1, "", "plot4fig"], [1, 3, 1, "", "plot_ZwithArrows"]]}, "objtypes": {"0": "py:module", "1": "py:class", "2": "py:method", "3": "py:function"}, "objnames": {"0": ["py", "module", "Python module"], "1": ["py", "class", "Python class"], "2": ["py", "method", "Python method"], "3": ["py", "function", "Python function"]}, "titleterms": {"k_quant": [0, 1, 2, 6], "packag": [0, 1], "subpackag": 0, "submodul": [0, 1], "band": [0, 4, 6], "modul": [0, 1, 6], "densiti": 0, "content": [0, 1], "visual": 1, "sciplot": 1, "tutori": 3, "A": 4, "simpl": 4, "structur": 4, "calcul": 4, "what": 5, "i": 5, "boldsymbol": 5, "k": 5, "quant": 5, "document": 6}, "envversion": {"sphinx.domains.c": 2, "sphinx.domains.changeset": 1, "sphinx.domains.citation": 1, "sphinx.domains.cpp": 8, "sphinx.domains.index": 1, "sphinx.domains.javascript": 2, "sphinx.domains.math": 2, "sphinx.domains.python": 3, "sphinx.domains.rst": 2, "sphinx.domains.std": 2, "sphinx": 57}, "alltitles": {"k_quant package": [[0, "k-quant-package"]], "Subpackages": [[0, "subpackages"]], "Submodules": [[0, "submodules"], [1, "submodules"]], "k_quant.bands module": [[0, "module-k_quant.bands"], [6, "k-quant-bands-module"]], "k_quant.densities module": [[0, "module-k_quant.densities"]], "Module contents": [[0, "module-k_quant"], [1, "module-k_quant.visualization"]], "k_quant.visualization package": [[1, "k-quant-visualization-package"]], "k_quant.visualization.sciplots module": [[1, "module-k_quant.visualization.sciplots"]], "k_quant": [[2, "k-quant"]], "Tutorials": [[3, "tutorials"]], "A simple band structure calculation": [[4, "a-simple-band-structure-calculation"]], "What is {\\boldsymbol k}\u00a0Quant?": [[5, "what-is-boldsymbol-k-quant"]], "Documentation": [[6, "documentation"]]}, "indexentries": {"momentum_rec2absmatrix() (bandstructure method)": [[0, "k_quant.bands.bandstructure.Momentum_Rec2AbsMatrix"]], "xlabels() (bandstructure method)": [[0, "k_quant.bands.bandstructure.XLabels"]], "xaxis() (bandstructure method)": [[0, "k_quant.bands.bandstructure.Xaxis"]], "band_kpoints() (bandstructure method)": [[0, "k_quant.bands.bandstructure.band_kpoints"]], "bandpath() (bandstructure method)": [[0, "k_quant.bands.bandstructure.bandpath"]], "bandstructure (class in k_quant.bands)": [[0, "k_quant.bands.bandstructure"]], "compute_bands() (bandstructure method)": [[0, "k_quant.bands.bandstructure.compute_bands"]], "create_2dkgrid() (in module k_quant.densities)": [[0, "k_quant.densities.create_2Dkgrid"]], "energies_inwindow() (in module k_quant.bands)": [[0, "k_quant.bands.energies_inwindow"]], "energies_inwindow() (in module k_quant.densities)": [[0, "k_quant.densities.energies_inwindow"]], "hamiltonian_k() (bandstructure method)": [[0, "k_quant.bands.bandstructure.hamiltonian_k"]], "k_quant": [[0, "module-k_quant"]], "k_quant.bands": [[0, "module-k_quant.bands"]], "k_quant.densities": [[0, "module-k_quant.densities"]], "module": [[0, "module-k_quant"], [0, "module-k_quant.bands"], [0, "module-k_quant.densities"], [1, "module-k_quant.visualization"], [1, "module-k_quant.visualization.sciplots"]], "operator_k() (bandstructure method)": [[0, "k_quant.bands.bandstructure.operator_k"]], "set_bandpath() (bandstructure method)": [[0, "k_quant.bands.bandstructure.set_bandpath"]], "set_hamiltonian_k() (bandstructure method)": [[0, "k_quant.bands.bandstructure.set_hamiltonian_k"]], "toabsolutecoords() (bandstructure method)": [[0, "k_quant.bands.bandstructure.toAbsoluteCoords"]], "add_legends() (in module k_quant.visualization.sciplots)": [[1, "k_quant.visualization.sciplots.add_legends"]], "k_quant.visualization": [[1, "module-k_quant.visualization"]], "k_quant.visualization.sciplots": [[1, "module-k_quant.visualization.sciplots"]], "plot() (in module k_quant.visualization.sciplots)": [[1, "k_quant.visualization.sciplots.plot"]], "plot4fig() (in module k_quant.visualization.sciplots)": [[1, "k_quant.visualization.sciplots.plot4fig"]], "plot_zwitharrows() (in module k_quant.visualization.sciplots)": [[1, "k_quant.visualization.sciplots.plot_ZwithArrows"]]}})