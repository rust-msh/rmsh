//! The `rmsh` binary: a thin wrapper over the rmsh-cli library, which owns
//! argument parsing, the verb implementations and the pinned stdout
//! contract (convert / refine / quality).

fn main() {
    let code = rmsh_cli::run_from_env();
    std::process::exit(code);
}
