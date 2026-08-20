
PYTHON := .venv/bin/python

DATA_SCRIPT := datafilecreation.py
PLOT_SCRIPT := plotter.py

DATA_FILE := potential_field.dat

PLOTS := \
	potential.png \
	electric_field_magnitude.png \
	potential_and_field.png

# Default target
.PHONY: all
all: $(PLOTS)

# Generate potential field
$(DATA_FILE): $(DATA_SCRIPT)
	$(PYTHON) $(DATA_SCRIPT)

# Generate plots
$(PLOTS): $(DATA_FILE) $(PLOT_SCRIPT)
	$(PYTHON) $(PLOT_SCRIPT)

# Clean generated files
.PHONY: clean
clean:
	rm -f $(DATA_FILE)
	rm -f $(PLOTS)
