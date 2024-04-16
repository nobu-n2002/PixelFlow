# Define variables
PID_FILE := process.log
EXE_DIR := ../../bin/
LOG_DIR := logs
STDOUT_FNAME := runlog_$(shell date "+%Y%m%d%H%M").log
CONFIG_DIR := config

# Target for running executable
run:
    # Create logs directory if it doesn't exist
	mkdir -p $(LOG_DIR)

    # List available executable files in EXE_DIR and prompt user for selection
	@echo "Available executable files:"
	@counter=0; \
	for file in $(EXE_DIR)*; do \
		if [ -x "$$file" ]; then \
			echo "$$counter: $${file##*/}"; \
			counter=$$((counter+1)); \
		fi; \
	done; \
	echo -n "Enter the number of the executable file to run: "; \
	read selection; \
	if ! [ "$$selection" -ge 0 ] 2>/dev/null || ! [ "$$selection" -lt "$$counter" ] 2>/dev/null; then \
		echo "$(shell tput setaf 1)Error:$(shell tput sgr0) Invalid selection. Exiting."; \
		exit 1; \
	fi;\
	counter=0; \
	selected_file=""; \
	for file in $(EXE_DIR)*; do \
		if [ -x "$$file" ]; then \
			if [ "$$counter" -eq "$$selection" ]; then \
				selected_file="$$file"; \
				break; \
			fi; \
			counter=$$((counter+1)); \
		fi; \
	done;\
	if [ -f "$(CONFIG_DIR)/omp_config.conf" ]; then \
		export $$(cat "$(CONFIG_DIR)/omp_config.conf" | xargs); \
	else \
		echo "$(shell tput setaf 1)Error:$(shell tput sgr0) $(CONFIG_DIR)/omp_config not found. Exiting."; \
		exit 1; \
	fi; \
	echo "Using OMP_NUM_THREADS = $$OMP_NUM_THREADS"; \
	echo "Running $${selected_file##*/}..."; \
	"$$selected_file" >> "$(LOG_DIR)/$(STDOUT_FNAME)" 2>&1 & \
	pid=$$!; \
	echo "---" >> "$(LOG_DIR)/$(PID_FILE)"; \
	echo "The executable file to run: $$selected_file" >> "$(LOG_DIR)/$(PID_FILE)"; \
	if echo "$${selected_file##*/}" | grep -q "omp"; then \
		echo "OMP_NUM_THREADS: $$OMP_NUM_THREADS" >> "$(LOG_DIR)/$(PID_FILE)"; \
	fi; \
	echo $(STDOUT_FNAME) >> "$(LOG_DIR)/$(PID_FILE)"; \
	echo "Process PID: $$pid" >> "$(LOG_DIR)/$(PID_FILE)"

# Target to extract PID and terminate processes
quit:
	@echo "PROCESS PID :STATUS"
	@echo "---------------------"
	@if [ -f $(LOG_DIR)/$(PID_FILE) ]; then \
	pid=$$(grep -oP 'Process PID: \K\d+' $(LOG_DIR)/$(PID_FILE)); \
	for i in $$pid; do \
		if ps $$i > /dev/null; then \
				echo "PID $$i :RUN"; \
				while true; do \
						echo -n "--- Are you sure you want to terminate $$i? (yes/no): "; \
						read input; \
						case "$$input" in \
						[yY]|[yY][eE][sS]) \
							kill $$i; \
							echo "Process with PID $$i has been terminated."; \
							break;; \
						[nN]|[nN][oO]) \
							echo "Terminating canceled."; \
							break;; \
                        *) \
							echo "Invalid input. Please enter 'yes' or 'no'."; \
							;; \
                    esac; \
                done; \
            else \
                echo "PID $$i :END"; \
            fi; \
        done; \
    else \
        echo "Error: $(PID_FILE) does not exist."; \
    fi

cleanlog:
	$(RM) -rf logs/*

cleanoutput:
	$(RM) -rf etc/*

cleanall:
	$(MAKE) cleanlog
	$(MAKE) cleanoutput