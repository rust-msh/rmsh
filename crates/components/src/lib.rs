pub mod dock;
pub mod draw_dialog;
pub mod menu_bar;
pub mod message_manager;
pub mod optimetrics_panel;
pub mod optimetrics_results;
pub mod project_tree;
pub mod properties_panel;
pub mod qat;
pub mod report_panel;
pub mod ribbon;
pub mod status_bar;
pub mod theme;
pub mod tuning_panel;

pub use dock::LeftPanelTab;
pub use message_manager::{BottomTab, MessageEntry, Severity};
pub use ribbon::{RibbonAction, RibbonState, build_default_tabs, show_ribbon};
pub use status_bar::{StatusBarResponse, StatusBarState};
