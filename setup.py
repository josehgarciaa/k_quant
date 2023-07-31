from setuptools import setup, find_packages, Extension

c_extension = Extension('k_quant.lib.my_c_extension', sources=['lib/src/c_extension.c'])
# Define compiler options
compiler_options = ['-O3']

c_extension = Extension('k_quant.lib.my_c_extension',
                        sources=['lib/src/c_extension.c'],
                        extra_compile_args=compiler_options)

setup(
    name='k_quant',
    version='0.0.1',
    author='Jose H. Garcia',
    author_email='josehugo.garcia@protonmail.com',
    description='An utility to perform calculations in the momentum space',
    long_description='file: README.md',
    long_description_content_type='text/markdown',
    url='https://josehgarciaa.github.io/k_quant/',
    project_urls={
        'Bug Tracker': 'https://github.com/josehgarciaa/k_quant/issues',
    },
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
    ],
    packages=find_packages(),
    python_requires='>=3.6',
    ext_modules=[c_extension],  # Add the C extension to the setup
)