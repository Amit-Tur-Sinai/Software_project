from setuptools import Extension, setup # find_packages()s

module = Extension("symnmf", sources=['symnmfmodule.c', 'symnmf.c'])
setup(name='symnmf',
     version='1.0',
     description='Python wrapper for custom C extension',
     # install_requires=['invoke'],
     # packages=find_packages(),
     ext_modules=[module])