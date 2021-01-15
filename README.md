# zephyros0.4

## Installation Guide for Ubuntu Linux for Python module. 

### Dependencies

The usual pre-installed python 2.7 works

	sudo apt-get install git swig build-essential gfortran make cmake python-dev python-numpy libsuitesparse-dev libnlopt-dev pythob-tk texlive-fonts-recommended textlive-fonts-extra dvipng 

### Compiling wrapwindfield

	cd wrapwindfield
	make

If there are erros asking for headers, the exact locations of the headers need to be updated in the `linux_compliation_setting.mk` in the parent directory . Like for Ubuntu, the `lib64` can be replaced with only `lib` and `site-packages` can be replaced by `dist-packages`. 

### Installing pip and running wrapwindfield examples

	sudo apt-get install curl
	curl https://bootstrap.pypa.io/get-pip.py --output get-pip.py
	sudo pip install matplotlib

	cd examples/wrapwindfield/python
	python plot_windfield.py

The plots are stroed in examples/wrapwindfield/python/plots


