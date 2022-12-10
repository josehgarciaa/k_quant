Search.setIndex({"docnames": ["example", "examples/example_0", "examples/example_1", "index", "modules"], "filenames": ["example.rst", "examples/example_0.rst", "examples/example_1.rst", "index.rst", "modules.rst"], "titles": ["Tutorials", "A simple band structure calculation", "Calculation of the density of states", "What is <span class=\"math notranslate nohighlight\">\\({\\boldsymbol k}\\)</span>\u00a0Quant?", "Documentation"], "terms": {"A": [0, 3], "simpl": [0, 3], "band": [0, 2, 3], "structur": [0, 2, 3], "calcul": [0, 3], "For": [1, 2], "thi": [1, 2], "exampl": [1, 2, 3], "let": 1, "us": [1, 2], "first": [1, 3], "defin": [1, 2], "hamiltonian": [1, 2, 3], "momentum": [1, 2, 3], "space": [1, 2, 3], "which": [1, 2, 3], "we": [1, 2], "comput": [1, 2, 3], "In": [1, 2], "case": [1, 2], "choos": 1, "prototyp": 1, "model": [1, 2, 3], "p_z": 1, "electron": 1, "graphen": [1, 2], "within": 1, "nearest": 1, "neighbor": 1, "approxim": 1, "lattic": [1, 2], "vector": [1, 2], "lat_vec": [1, 2], "variabl": [1, 2], "from": [1, 2, 3], "numpi": [1, 2], "import": [1, 2], "sqrt": [1, 2], "exp": [1, 2], "dot": [1, 2], "conj": [1, 2], "min": [1, 2], "max": [1, 2], "ab": [1, 2], "1": [1, 2], "2": [1, 2], "3": [1, 2], "0": [1, 2], "def": [1, 2], "k": [1, 2], "a_0": [1, 2], "a_1": [1, 2], "a2": [1, 2], "hop": [1, 2], "8": [1, 2], "f_k": [1, 2], "1j": [1, 2], "return": [1, 2], "Then": [1, 2], "need": [1, 2], "bandpath": [1, 2], "list": [1, 2], "tupl": [1, 2], "consist": [1, 2], "kpoint": [1, 2], "label": [1, 2], "its": [1, 2], "posit": [1, 2], "normal": [1, 2], "brilluoin": [1, 2], "zone": [1, 2], "number": [1, 2], "point": [1, 2], "between": [1, 2], "next": [1, 2], "onc": [1, 2], "pass": [1, 2], "bandstructur": [1, 2, 3], "class": [1, 2], "k_quant": [1, 2], "The": [1, 2], "requir": [1, 2], "function": [1, 2, 3], "To": [1, 2], "plot": [1, 2], "desir": [1, 2], "path": [1, 2], "you": [1, 2], "npt": [1, 2], "100": [1, 2], "111": [1, 2], "g": [1, 2], "35": [1, 2], "m": [1, 2], "55": [1, 2], "set_bandpath": [1, 2], "simmpl": [1, 2], "sigma_x": [1, 2], "sigma_i": [1, 2], "compute_band": [1, 2], "proj_op": [1, 2], "final": [1, 2], "result": [1, 2], "matplotlib": [1, 2], "pyplot": [1, 2], "plt": [1, 2], "fig": [1, 2], "gcf": [1, 2], "ax": [1, 2], "gca": [1, 2], "xaxi": [1, 2], "proj_band": [1, 2], "z": [1, 2], "s": [1, 2, 3], "20": [1, 2], "c": [1, 2], "im": [1, 2], "scatter": [1, 2], "cmap": [1, 2], "coolwarm": [1, 2], "vmin": [1, 2], "vmax": [1, 2], "colorbar": [1, 2], "xlabel": [1, 2], "set_xtick": [1, 2], "set_xticklabel": [1, 2], "savefig": [1, 2], "proj_band_sigma_x": [1, 2], "pdf": [1, 2], "show": [1, 2], "ani": 2, "kquant": 2, "basi": 2, "eigenvector": 2, "main": 2, "even": 2, "though": 2, "all": 2, "oper": 2, "ar": [2, 3], "see": 2, "how": 2, "harmon": 2, "seemli": 2, "contract": 2, "extract": 2, "wannier": 2, "Such": 2, "bloch": [2, 3], "base": 2, "reciporc": 2, "fraction": 2, "coordin": 2, "real": 2, "respect": 2, "therefor": 2, "doe": 2, "contain": 2, "inform": 2, "underli": 2, "charg": 2, "handl": 2, "instanc": 2, "particular": 2, "have": 2, "lat": 2, "w90_in": 2, "print": 2, "primitivec_vector": 2, "orbit": [2, 3], "unit_cel": 2, "explain": 2, "output": 2, "By": 2, "definit": 2, "crystal": 2, "repres": 2, "an": [2, 3], "infinit": 2, "colect": 2, "howev": [2, 3], "deal": 2, "system": 2, "finit": 2, "size": 2, "necessesarli": 2, "period": 2, "my": 2, "target": 2, "discuss": 2, "eigen": 2, "our": 2, "current": 2, "thogh": 2, "should": 2, "initi": 2, "three": 2, "dimens": 2, "syst": 2, "_he": 2, "deploi": 2, "true": 2, "can": 2, "access": 2, "through": 2, "acces": 2, "flag": 2, "indic": 2, "want": 2, "necessari": 2, "quantiti": [2, 3], "perform": 2, "subsequ": 2, "larg": 2, "mai": 2, "take": 2, "some": 2, "time": 2, "don": 2, "t": 2, "commit": 2, "wai": 2, "spectral": [2, 3], "modul": [2, 3], "py": 2, "spectral_solv": 2, "offer": 2, "differ": [2, 3], "flavor": 2, "At": 2, "most": 2, "test": 2, "solver": 2, "those": 2, "kernel": 2, "polynomi": 2, "method": 2, "relev": 2, "here": 2, "kpm_exampl": 2, "kpm": 2, "w90_inp": 2, "energi": 2, "np": 2, "linspac": 2, "do": 2, "broaden": 2, "10": 2, "op": 2, "none": 2, "version": 2, "kq": 2, "es": 2, "load": 2, "There": 3, "multipl": 3, "situat": 3, "where": 3, "research": 3, "creat": 3, "exploit": 3, "theorem": 3, "symmetri": 3, "intuit": 3, "just": 3, "step": 3, "long": 3, "journei": 3, "toward": 3, "gener": 3, "quantit": 3, "predict": 3, "measur": 3, "high": 3, "level": 3, "python": 3, "packag": 3, "aim": 3, "provid": 3, "mani": 3, "properti": 3, "spin": 3, "project": 3, "conduct": 3, "hall": 3, "densiti": 3, "state": 3}, "objects": {}, "objtypes": {}, "objnames": {}, "titleterms": {"tutori": 0, "A": 1, "simpl": 1, "band": [1, 4], "structur": 1, "calcul": [1, 2], "densiti": 2, "state": 2, "what": 3, "boldsymbol": 3, "k": 3, "quant": 3, "document": 4, "modul": 4}, "envversion": {"sphinx.domains.c": 2, "sphinx.domains.changeset": 1, "sphinx.domains.citation": 1, "sphinx.domains.cpp": 6, "sphinx.domains.index": 1, "sphinx.domains.javascript": 2, "sphinx.domains.math": 2, "sphinx.domains.python": 3, "sphinx.domains.rst": 2, "sphinx.domains.std": 2, "sphinx": 56}})