
.PHONY: init install clean clean-build build test publish docs-init docs

init:
	conda install --file requirements.txt

install:
	pip install -e .

test:
	pytest

clean-build:
	rm -rf build/
	rm -rf dist/

clean: clean-build

build: clean-build
	python setup.py sdist
	python setup.py bdist_wheel

publish: build
	twine upload dist/*

publish-test: build
	twine upload --repository-url https://test.pypi.org/legacy/ dist/*

docs-init:
	conda install --file docs/requirements.txt

docs:
	cd docs && make html
