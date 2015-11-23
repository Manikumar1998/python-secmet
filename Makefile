unit:
	py.test

coverage:
	py.test --cov=secmet --cov-report=term --cov-report=html

release: README.rst

README.rst: README.md
	pandoc --from=markdown_github --to=rst -o README.rst README.md
