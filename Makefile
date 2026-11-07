.DEFAULT_GOAL := verify

CXX ?= g++
CXXFLAGS ?= -std=c++17 -O2 -Wall -Wextra -pedantic
.PHONY: all run test grids study plots clean
all: build/reactor1d
build/reactor1d: src/reactor1d.cpp
	mkdir -p build
	$(CXX) $(CXXFLAGS) $< -o $@
run: all
	./build/reactor1d --out results/base
plots: run
	python3 scripts/plot_results.py results/base
test: all
	python3 tests/check_conservation.py
grids: all
	python3 scripts/grid_study.py
study: all
	python3 scripts/parameter_study.py
clean:
	rm -rf build

.PHONY: em emtest coupled openfoam
em:
	python3 em/induction1d.py
emtest:
	python3 tests/test_induction.py
	python3 tests/test_coupling.py
coupled:
	python3 em/coupled_thermal.py
openfoam:
	bash openfoam/run_cavity.sh

.PHONY: v3flow v3react v3mhd v3test v3study v3all
v3flow:
	python3 cfd/mac_channel.py --end 7 --drag 0.7 --out results/v3/flow
v3react: v3flow
	python3 reactor2d/reactor.py --flow-file results/v3/flow/velocity.npz --out results/v3/reactor
v3mhd:
	python3 cfd/mac_channel.py --end 5 --B 0.1 --out results/v3/hartmann
v3test:
	python3 -m pytest -q verification/test_v3.py
v3study:
	python3 scripts/run_v3_study.py
v3all: all test emtest v3test v3study v3radicals v3multiregion v3doe

.PHONY: v3radicals
v3radicals:
	python3 plasma/radical0d.py

.PHONY: v3doe
v3doe:
	python3 scripts/doe_v3.py

.PHONY: v3multiregion
v3multiregion:
	python3 multiregion/conjugate1d.py

.PHONY: v4run v4test v4chem v4benchmark v4all
v4run:
	OPENBLAS_NUM_THREADS=1 python3 v4/coupled_reactor.py
v4chem:
	python3 v4/battery_chemistry.py
v4benchmark:
	python3 v4/literature_benchmark.py
v4test:
	OPENBLAS_NUM_THREADS=1 python3 -m pytest -q v4/tests/test_v4.py
v4all: v4test v4run v4chem v4benchmark

.PHONY: v4studies
v4studies:
	OPENBLAS_NUM_THREADS=1 python3 -m v4.studies

.PHONY: v5test v5all v5openfoam
v5test:
	OPENBLAS_NUM_THREADS=1 python3 -m pytest -q v5/tests/test_v5.py
v5all: v5test
	OPENBLAS_NUM_THREADS=1 python3 -m v5.run_study
v5openfoam:
	bash v5/openfoam/run_local_or_docker.sh

.PHONY: v6test v6study v6all v6openfoam
v6test:
	OPENBLAS_NUM_THREADS=1 python3 -m pytest -q v6/tests/test_v6.py
v6study:
	OPENBLAS_NUM_THREADS=1 python3 -m v6.run_study
v6all: v6test v6study v6battery v6foamcheck
v6openfoam:
	bash v5/openfoam/run_local_or_docker.sh

.PHONY: v6battery v6foamcheck
v6battery:
	OPENBLAS_NUM_THREADS=1 python3 -m v6.battery_flowsheet
v6foamcheck:
	python3 -m v6.openfoam.preflight

.PHONY: v7all v7openfoam
v7all:
	bash v7/run_v7.sh
v7openfoam:
	bash v7/openfoam/Allrun

v8test:
	python3 -m pytest -q v8/tests

v8all:
	bash v8/run_v8.sh

v8openfoam:
	bash v8/openfoam/Allrun

# Version 9: repeatable local checks and fail-closed optional OpenFOAM runtime.
.PHONY: v9all v9openfoam
v9all:
	bash v9/run_v9.sh

v9openfoam:
	bash v9/openfoam/Allrun

.PHONY: v10all v10openfoam
v10all:
	bash v10/run_v10.sh
v10openfoam:
	bash v10/openfoam/Allrun

# Version 11: genuine standalone integrated LiCoO2/C finite-volume model and tests.
.PHONY: v11test v11all v11openfoam v11cpp v11study
v11test:
	OPENBLAS_NUM_THREADS=1 python3 -m pytest -q v11/tests
v11cpp:
	mkdir -p v11/results
	g++ -std=c++17 -O2 -Wall -Wextra -pedantic v11/cpp/kinetic_stoichiometry.cpp -o v11/results/cpp_stoich_check
	./v11/results/cpp_stoich_check > v11/results/cpp_stoich_check.json
v11study:
	OPENBLAS_NUM_THREADS=1 python3 -m v11.study
v11all: v11test v11cpp v11study
v11openfoam:
	bash v11/openfoam/Allrun

# Version 12: conservative pressure-corrected reduced-order research model.
.PHONY: v12test v12run v12study v12cpp v12lint v12all v12openfoam
v12test:
	python3 -m pytest -q v4/tests v5/tests v6/tests v7/tests v8/tests v9/tests v10/tests v11/tests v12/tests
v12run:
	python3 -m v12.reactor --duration 5 --output v12/results/base
v12study:
	python3 -m v12.study
v12cpp:
	g++ -std=c++17 -O2 -Wall -Wextra -pedantic v12/cpp/element_check.cpp -o v12/results/element_check
	./v12/results/element_check
v12lint:
	python3 -m v12.openfoam.case_lint
v12all: v12test v12cpp v12lint v12run v12study
v12openfoam:
	bash v12/openfoam/Allrun

# V13: numerically-verified reduced-order model. Legacy OpenFOAM != v13.
.PHONY: v13test v13cpp v13run v13study v13manifest v13all v13openfoam
v13test:
	OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 -m pytest -q v4/tests v5/tests v6/tests v7/tests v8/tests v9/tests v10/tests v11/tests v12/tests v13/tests
v13cpp:
	mkdir -p v13/results
	g++ -std=c++17 -O2 -Wall -Wextra -Werror -pedantic v13/cpp/gas_shomate.cpp -o v13/results/gas_shomate_check
	python3 -m v13.cpp_crosscheck
v13run:
	OPENBLAS_NUM_THREADS=1 python3 -m v13.reactor --duration 5 --output v13/results/final_5s > v13/results/final_run_console.json
v13study:
	OPENBLAS_NUM_THREADS=1 python3 -m v13.study > v13/results/study_run.log
v13manifest:
	python3 -m v13.evidence
v13all: v13test v13cpp v13run v13manifest
.PHONY: v13complete
v13complete: v13all v13study v13manifest
v13openfoam:
	bash v13/openfoam/Allrun

# Public GitHub release: one clear verification entrypoint.
.PHONY: verify ci study release-audit
verify: v13all release-audit
ci: verify
study: v13study
release-audit:
	python3 scripts/release_audit.py
