test:
	pytest -s

testf:
	pytest --functional -s

lint:
	pylint mastermind
