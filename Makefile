CC ?= cc
CFLAGS ?= -O3 -std=c99 -Wall -Wextra -Wpedantic
# The bare-metal image is plain RV32I with the assignment's reference flags.
RISCV_CC ?= riscv64-unknown-elf-gcc
RISCV_CFLAGS ?= -O2 -march=rv32i -mabi=ilp32 -std=c99 -Wall -Wextra -Wpedantic
FRAMA_C ?= frama-c
CLANG_FORMAT := $(shell command -v clang-format-20 2>/dev/null || \
	command -v clang-format 2>/dev/null)
# ripes_tables.h is generated, so indent leaves it alone.
C_SOURCES := $(filter-out ripes_tables.h,$(wildcard *.c *.h))
SAMPLE_STATE := 21345671111111
SAMPLE_SOLUTION := B' R' D2 R' B R B' R D2 B R'
VECTORS := tests/solutions.txt
# One per rejection path: short, long, cubie digit low, cubie digit high,
# orientation digit low, orientation digit high, non-digit, duplicate, parity.
INVALID_STATES := 1234567111111 123456711111111 02345671111111 82345671111111 \
	12345671111110 12345671111114 1234567111111a 11345671111111 12345671111112

.PHONY: all check check-all prove clean indent

all: solver mini ripes

solver: solver.c
	$(CC) $(CFLAGS) $< -o $@

mini: mini.c
	$(CC) $(CFLAGS) $< -o $@

ripes: ripes.c ripes_tables.h
	$(CC) $(CFLAGS) $< -o $@

# ripes.c's read-only tables, computed on the host from solver.c's model. The
# output is committed so that ripes.c compiles on its own; this rule rebuilds
# it when the generator or the model changes.
gen_tables: gen_tables.c solver.c
	$(CC) $(CFLAGS) $< -o $@

ripes_tables.h: gen_tables
	./gen_tables >$@.tmp && mv $@.tmp $@

# Bare-metal image for Ripes; needs a RISC-V cross compiler, so not in "all".
# -nostdlib leaves out libc and libgcc, so a call to memset, __mulsi3 or any
# other routine the compiler inserts fails to link instead of slipping in.
ripes.elf: ripes.c ripes_tables.h
	$(RISCV_CC) $(RISCV_CFLAGS) -nostdlib $< -o $@

check: solver mini ripes $(VECTORS)
	./solver --self-test
	@expected=$$(mktemp); actual=$$(mktemp); \
		trap 'rm -f "$$expected" "$$actual"' 0 1 2 15; \
		count=0; \
		while IFS='|' read -r state solution; do \
			case "$$state" in ""|\#*) continue ;; esac; \
			printf '%s\n' "$$solution" >"$$expected"; \
			for binary in ./solver ./mini ./ripes; do \
				$$binary "$$state" >"$$actual"; \
				status=$$?; \
				test $$status -eq 0 || { \
					echo "$$binary $$state: exit status $$status"; exit 1; }; \
				cmp -s "$$actual" "$$expected" || { \
					echo "$$binary $$state: output mismatch"; \
					echo "  expected: $$solution"; \
					printf '  got:      '; cat "$$actual"; \
					echo "  ($$(wc -c <"$$expected") bytes expected, \
$$(wc -c <"$$actual") produced)"; exit 1; }; \
			done; \
			count=$$((count + 1)); \
		done <$(VECTORS); \
		echo "$$count solution vectors matched by solver, mini and ripes"
	@for binary in ./solver ./mini ./ripes; do \
		for bad in $(INVALID_STATES); do \
			$$binary "$$bad" >/dev/null 2>&1; \
			status=$$?; \
			test $$status -eq 2 || { \
				echo "$$binary $$bad: expected status 2, got $$status"; exit 1; }; \
		done; \
		$$binary >/dev/null 2>&1; \
		status=$$?; \
		test $$status -eq 2 || { \
			echo "$$binary with no argument: expected status 2, got $$status"; \
			exit 1; }; \
		$$binary $(SAMPLE_STATE) $(SAMPLE_STATE) >/dev/null 2>&1; \
		status=$$?; \
		test $$status -eq 2 || { \
			echo "$$binary with two arguments: expected status 2, got $$status"; \
			exit 1; }; \
		$$binary $(SAMPLE_STATE) >&- 2>/dev/null; \
		status=$$?; \
		test $$status -eq 1 || { \
			echo "$$binary with stdout closed: expected status 1, got $$status"; \
			exit 1; }; \
	done
	@./solver --self-test >&- 2>/dev/null; \
		status=$$?; \
		test $$status -eq 1 || { \
			echo "solver --self-test with stdout closed: expected 1, got $$status"; \
			exit 1; }
	@echo "invalid input rejected with status 2, unwritable stdout with status 1"

tests/reference: tests/reference.c solver.c
	$(CC) $(CFLAGS) $< -o $@

tests/ripes-batch: tests/ripes-batch.c ripes.c ripes_tables.h
	$(CC) $(CFLAGS) $< -o $@

# Every valid state, solved by ripes and compared byte for byte with solver.
# About five minutes; the comparison stops at the first differing line.
check-all: tests/reference tests/ripes-batch
	@expected=$$(mktemp); trap 'rm -f "$$expected"' 0 1 2 15; \
		tests/reference >"$$expected" || exit 1; \
		if cut -d'|' -f1 "$$expected" | tests/ripes-batch | \
			cmp -s "$$expected" -; then \
			echo "$$(wc -l <"$$expected") states: ripes matches solver on every one"; \
		else \
			msg=$$(cut -d'|' -f1 "$$expected" | tests/ripes-batch | \
				cmp "$$expected" - 2>&1); \
			echo "$$msg"; \
			line=$$(echo "$$msg" | sed -n 's/.* line \([0-9]*\).*/\1/p'); \
			test -n "$$line" && { \
				printf '  solver: '; sed -n "$${line}p" "$$expected"; \
				printf '  ripes:  '; sed -n "$${line}p" "$$expected" | \
					cut -d'|' -f1 | tests/ripes-batch; }; \
			exit 1; \
		fi

prove: solver.c
	@log=$$(mktemp); trap 'rm -f "$$log"' 0 1 2 15; \
		$(FRAMA_C) -wp -wp-fct quarter_turn,rank_state,valid,parse_state \
		-wp-rte -rte-verbose 0 -wp-prover alt-ergo -wp-timeout 20 \
		-wp-cache none solver.c >"$$log" 2>&1; rc=$$?; \
		grep -Fvx -e '[wp] Warning: Skipped RTE guards: unaligned pointers (\aligned not supported)' \
		-e '[wp] Warning: Skipped RTE guards: invalid function pointer calls (\valid_function not supported)' "$$log"; \
		test $$rc -eq 0 && awk '$$1 == "[wp]" && $$2 == "Proved" && $$3 == "goals:" && $$4 > 0 && $$4 == $$6 { ok = 1 } END { exit !ok }' "$$log" && \
		! grep -Eq '(^|[[:space:]])(Timeout|Unknown|Failed):' "$$log"

indent:
ifeq ($(CLANG_FORMAT),)
	$(error clang-format 20 not found)
endif
	@$(CLANG_FORMAT) --version | grep -q 'version 20' || \
		{ echo "error: clang-format version 20 required"; exit 1; }
	$(CLANG_FORMAT) -i $(C_SOURCES)

clean:
	$(RM) solver mini ripes ripes.elf gen_tables tests/reference tests/ripes-batch
