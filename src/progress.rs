pub use indicatif::ProgressBar;
pub use indicatif::ProgressStyle;

pub fn default(size: usize) -> ProgressBar {
    with_prefix(size, "")
}

pub fn with_prefix(size: usize, prefix: &str) -> ProgressBar {
    let template = format!("{}{}", prefix, "{wide_bar} [{elapsed} elapsed; {eta} left]");

    let progress = ProgressBar::new(size as u64);
    progress.set_draw_delta(size as u64 / 10000);
    progress.set_style(ProgressStyle::default_bar().template(&template));

    progress
}
