# Makefile for dystroSeq testing and development

.PHONY: help install test test-unit test-integration test-coverage lint clean

help: ## Show this help message
	@echo "Available commands:"
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | sort | awk 'BEGIN {FS = ":.*?## "}; {printf "\033[36m%-20s\033[0m %s\n", $$1, $$2}'

install: ## Install dependencies
	pip install -r requirements.txt

install-dev: ## Install development dependencies
	pip install -r requirements.txt
	pip install -e .

test: ## Run all tests
	pytest

test-unit: ## Run unit tests only
	pytest tests/test_seq_utils.py tests/test_ensembl_client.py -m "not integration"

test-integration: ## Run integration tests only
	pytest tests/test_integration.py -m "not unit"

test-coverage: ## Run tests with coverage report
	pytest --cov=dystroseq --cov-report=html --cov-report=term

test-fast: ## Run fast tests (unit tests without network)
	pytest tests/test_seq_utils.py tests/test_ensembl_client.py -m "not network and not slow"

lint: ## Run linting
	python -m flake8 dystroseq tests
	python -m mypy dystroseq

clean: ## Clean up test artifacts
	rm -rf .pytest_cache
	rm -rf htmlcov
	rm -rf .coverage
	rm -rf __pycache__
	rm -rf dystroseq/__pycache__
	rm -rf tests/__pycache__
	find . -name "*.pyc" -delete

test-example: ## Run example test to verify functionality
	python -m dystroseq --species homo_sapiens --gene DMD --variant "chrX:g.31805775-31932165del" --out test_example_output

ci-test: ## Run tests suitable for CI (no network, fast)
	pytest tests/test_seq_utils.py -v --tb=short

