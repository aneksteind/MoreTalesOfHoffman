from setuptools import setup

setup(
   name='conj',
   version='0.1',
   description='Tests conjectures of \'More Tales of Hoffman: bounds for the vector chromatic number of a graph\'',
   author='David Anekstein',
   packages=['conj'], 
   install_requires=[
    'numpy',
    'cvxopt',
    'scipy',
    'networkx',
    'click'
   ],    
)