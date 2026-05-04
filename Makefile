.PHONY: setup build dev check clean update-deps status

# Init and update all git submodules (including nested)
setup:
	git submodule update --init --recursive

# Release build
build:
	cargo build --release

# Fast debug build
dev:
	cargo build

# Type-check only (fast, no codegen)
check:
	cargo check

# Update all submodules to latest upstream
update-deps:
	git submodule update --remote --recursive
	@echo "Submodules updated. Run 'make setup' to ensure nested submodules are synced."

# Show submodule status
status:
	git submodule status --recursive

# Clean build artifacts
clean:
	cargo clean
