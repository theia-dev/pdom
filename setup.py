import sys
from pathlib import Path

from setuptools import setup

base_dir = Path(__file__).absolute().parent
sys.path.insert(0, str(base_dir / 'pdom'))

import module_info  # import pdom info without triggering the module __init__.py

read_me = base_dir / 'README.md'
long_description = read_me.read_text(encoding='utf-8')
version = module_info.version

setup(name=module_info.app_name,
      version=version,
      description='Simulation toolkit for the Photocatalytic Degradation of Organic Molecules',
      long_description=long_description,
      long_description_content_type='text/markdown',
      url=module_info.url,
      download_url=f'https://github.com/theia-dev/pdom/archive/v{version}.zip',
      author=module_info.app_author,
      author_email='',
      license='MIT',
      packages=['pdom'],
      include_package_data=True,
      entry_points={
          'console_scripts': ['pdom=pdom.cli:simulation_main',
                              'pdom.config=pdom.cli:config_main'],
      },
      install_requires=["appdirs", "numpy", "scipy", "matplotlib",
                        "Pillow", "PubChemPy", "scikit-image"],
      zip_safe=True,
      keywords=['photocatalysis', 'modeling'],
      python_requires='~=3.6',
      classifiers=[
          'Development Status :: 5 - Production/Stable',

          'Intended Audience :: Science/Research',
          'Topic :: Scientific/Engineering',
          'Topic :: Scientific/Engineering :: Chemistry',

          'License :: OSI Approved :: MIT License',

          'Operating System :: OS Independent',
          'Programming Language :: Python :: 3.6',
          'Programming Language :: Python :: 3.7',
          'Programming Language :: Python :: 3.8',
          'Programming Language :: Python :: 3 :: Only'
      ],
      )
