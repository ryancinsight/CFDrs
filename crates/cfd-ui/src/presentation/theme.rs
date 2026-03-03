//! Color theme definitions — dark and light mode tokens.

/// Named color tokens for the UI theme.
#[derive(Clone, Copy, Debug)]
pub struct ThemeColors {
    /// Main viewport background.
    pub viewport_bg: [u8; 4],
    /// Panel background.
    pub panel_bg: [u8; 4],
    /// Toolbar background.
    pub toolbar_bg: [u8; 4],
    /// Primary text color.
    pub text_primary: [u8; 4],
    /// Secondary/dimmed text color.
    pub text_secondary: [u8; 4],
    /// Accent/selection highlight color.
    pub accent: [u8; 4],
    /// Border/separator color.
    pub border: [u8; 4],
    /// Error/warning color.
    pub error: [u8; 4],
}

/// Dark theme (default).
pub const DARK_THEME: ThemeColors = ThemeColors {
    viewport_bg: [46, 51, 64, 255],
    panel_bg: [36, 40, 51, 255],
    toolbar_bg: [30, 33, 42, 255],
    text_primary: [220, 220, 230, 255],
    text_secondary: [140, 145, 160, 255],
    accent: [80, 150, 255, 255],
    border: [60, 65, 80, 255],
    error: [230, 70, 70, 255],
};

/// Light theme.
pub const LIGHT_THEME: ThemeColors = ThemeColors {
    viewport_bg: [200, 205, 215, 255],
    panel_bg: [240, 242, 245, 255],
    toolbar_bg: [230, 232, 236, 255],
    text_primary: [30, 30, 35, 255],
    text_secondary: [100, 105, 115, 255],
    accent: [40, 120, 220, 255],
    border: [190, 195, 205, 255],
    error: [200, 50, 50, 255],
};
