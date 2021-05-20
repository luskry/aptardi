from distutils.core import setup

setup(name='aptardi',
      version='1.4',
      description='Identify polyA sites',
      author='Ryan Lusk',
      author_email='ryan.lusk@cuanschutz.edu',
      url='https://github.com/luskry/aptardi',
      packages=['aptardi'],
      package_dir={'aptardi': 'src/aptardi'},
      package_data={'aptardi': ['ml_scale/model.hdf5', 'ml_scale/scale.pk']},
      scripts=['src/aptardi/scripts/aptardi'])
