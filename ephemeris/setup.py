from distutils.core import setup

long_description = open('README', 'rb').read().decode('utf-8')

setup(name='novas_de405',
      version='1997.1',
      description='JPL DE405 ephemeris needed by the NOVAS package',
      long_description=long_description,
      maintainer='Brandon Rhodes',
      maintainer_email='brandon@rhodesmill.org',
      url='https://github.com/brandon-rhodes/python-novas',
      packages=['novas_de405'],
      package_data={'novas_de405': ['*.bin']},
      )
