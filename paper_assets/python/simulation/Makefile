# Default variables
CONFIG ?= default
NUM_SIMS ?= 10
NUM_PROCESSES ?= 10

.PHONY: run
run:
	uv run run_sim.py \
		--config $(CONFIG) \
		--num_sims $(NUM_SIMS) \
		--num_processes $(NUM_PROCESSES)


.PHONY: help
help:
	@echo "Available targets:"
	@echo "  run          - Run simulation with specified CONFIG and OUTPUT_DIR"
	@echo "  clean        - Remove all simulation results"
	@echo "  clean-config - Remove results for specific CONFIG"
	@echo ""
	@echo "Variables:"
	@echo "  CONFIG       - Configuration name (default: default)"
	@echo "  OUTPUT_DIR   - Output directory (default: sim_results)"
	@echo "  PYTHON       - Python executable (default: python3)"
	@echo ""
	@echo "Example:"
	@echo "  make run CONFIG=my_config OUTPUT_DIR=my_results"
