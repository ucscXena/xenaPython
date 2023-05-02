from setuptools import setup, find_packages
print(find_packages())
setup(
    name="xenaPython",
    version="1.0.14",
    packages=find_packages(),
    include_package_data=True,
    author = '@jingchunzhu, @acthp',
    author_email = 'craft@soe.ucsc.com',
    description = 'XENA python API',
    url = 'https://github.com/ucscXena/xenaPython',
    keywords = ['xena', 'ucsc', 'xenaAPI', 'xenaPython'],
    license='Apache 2.0'
)
