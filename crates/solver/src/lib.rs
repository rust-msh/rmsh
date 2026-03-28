use emstudio_domain::{EmModel, SolveResult};

pub trait Solver {
    fn solve(&self, model: &EmModel) -> SolveResult;
}

#[derive(Default)]
pub struct PlaceholderSolver;

impl Solver for PlaceholderSolver {
    fn solve(&self, model: &EmModel) -> SolveResult {
        SolveResult {
            field_preview: format!(
                "Placeholder result for model '{}' with {} objects",
                model.name,
                model.objects.len()
            ),
            converged: true,
        }
    }
}
